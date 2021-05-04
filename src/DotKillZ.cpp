#include "DotKill.h"

template <typename T>
void applyFieldBlend(PVideoFrame& srcc, PVideoFrame& srcn, PVideoFrame& dst, int order, int planecount)
{
    const int plane[3] = { PLANAR_Y, PLANAR_U, PLANAR_V };

    for (int i = 0; i < planecount; ++i)
    {
        const int stride = dst->GetPitch(plane[i]) / sizeof(T);
        const int width = srcc->GetRowSize(plane[i]) / sizeof(T);
        const int height = srcc->GetHeight(plane[i]);
        const T* srccp = reinterpret_cast<const T*>(srcc->GetReadPtr(plane[i]));
        const T* srcnp = reinterpret_cast<const T*>(srcn->GetReadPtr(plane[i]));
        T* dstp = reinterpret_cast<T*>(dst->GetWritePtr(plane[i]));

        if (order)
        {
            srccp += stride;
            srcnp += stride;
            dstp += stride;
        }

        for (int h = order; h < height; h += 2)
        {
            for (int w = 0; w < width; ++w)
                dstp[w] = (srccp[w] + srcnp[w] + 1) / 2;

            srccp += 2LL * stride;
            srcnp += 2LL * stride;
            dstp += 2LL * stride;
        }
    }
}

template <typename T>
void applyDotcrawInverse(PVideoFrame& srcc, PVideoFrame& srcn, PVideoFrame& dst, int order, int depth, int planecount)
{
    const int plane[3] = { PLANAR_Y, PLANAR_U, PLANAR_V };
    const int shift = depth - 8;
    const int lmax = 235 << shift;
    const int lmin = 16 << shift;

    for (int i = 0; i < planecount; ++i)
    {
        const int stride = dst->GetPitch(plane[i]) / sizeof(T);
        const int width = srcc->GetRowSize(plane[i]) / sizeof(T);
        const int height = srcc->GetHeight(plane[i]);
        const T* srccp = reinterpret_cast<const T*>(srcc->GetReadPtr(plane[i]));
        const T* srcnp = reinterpret_cast<const T*>(srcn->GetReadPtr(plane[i]));
        T* dstp = reinterpret_cast<T*>(dst->GetWritePtr(plane[i]));

        if (order)
        {
            srccp += stride;
            srcnp += stride;
            dstp += stride;
        }

        for (int h = order; h < height; h += 2)
        {
            for (int w = 0; w < width; ++w)
            {
                dstp[w] = (srccp[w] + srcnp[w] + 1) / 2;

                if (h > 1)
                {
                    T l0val = dstp[w - 2 * stride];
                    T l2val = dstp[w];
                    int l0diff = dstp[w - 2 * stride] - srccp[w - 2 * stride];
                    int l2diff = dstp[w] - srccp[w];
                    if (i == 0)
                        dstp[w - stride] = std::clamp(srccp[w - stride] + (order ? l0diff : l2diff), lmin, lmax);
                    else
                        dstp[w - stride] = (l0val + l2val + 1) / 2; // simply use some kind of interpolation and discard one field?
                }
            }

            srccp += 2LL * stride;
            srcnp += 2LL * stride;
            dstp += 2LL * stride;
        }
    }
}

DotKillZ::DotKillZ(PClip _child, int _order, int _offset, IScriptEnvironment* env)
    : GenericVideoFilter(_child), order(_order), offset(_offset)
{
    if (vi.BitsPerComponent() == 32 || vi.IsRGB() || !vi.IsPlanar())
        env->ThrowError("DotKillZ: only YUV 8..16-bit planar format supported.");
    if (order < 0 || order > 1)
        env->ThrowError("DotKillZ: order must be either 0 or 1.");
    if (offset < 0 || offset > 4)
        env->ThrowError("DotKillZ: offset must be between 0..4.");

    if (vi.ComponentSize() == 1)
    {
        applyfb = applyFieldBlend<uint8_t>;
        applydi = applyDotcrawInverse<uint8_t>;
    }
    else
    {
        applyfb = applyFieldBlend<uint16_t>;
        applydi = applyDotcrawInverse<uint16_t>;
    }
}

PVideoFrame __stdcall DotKillZ::GetFrame(int n, IScriptEnvironment* env)
{
    PVideoFrame srcp = child->GetFrame(std::max(n - 1, 0), env);
    PVideoFrame srcc = child->GetFrame(n, env);
    PVideoFrame srcn = child->GetFrame(std::min(n + 1, vi.num_frames - 1), env);
    PVideoFrame dst = srcc;
    env->MakeWritable(&dst);

    const int planecount = std::min(vi.NumComponents(), 3);

    /*
    FIELD OFFSETS
    -1  0  1  2  3
    A1 B1 B1 C1 D1
    A2 B2 C2 D2 D2
    */

    if ((n + offset) % 5 == 0)
    {
        // current and next field are duplicates, complement field is from the same frame so do dotcrawl inverse on that as well
        applydi(srcc, srcn, dst, order, vi.BitsPerComponent(), planecount);
    }
    else if ((n + offset) % 5 == 1)
    {
        // current and previous field are duplicates so blend them together
        applyfb(srcc, srcp, dst, order, planecount);
    }
    else if ((n + offset) % 5 == 2)
    {
        // current and next complement field are duplicates so blend them together
        applyfb(srcc, srcn, dst, !order, planecount);
    }
    else if ((n + offset) % 5 == 3)
    {
        // current and previous field are duplicates, complement field is from the same frame so do dotcrawl inverse on that as well
        applydi(srcc, srcp, dst, !order, vi.BitsPerComponent(), planecount);
    }

    return dst;
}
