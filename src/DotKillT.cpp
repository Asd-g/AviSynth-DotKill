#include <vector>

#include "DotKill.h"

constexpr int blockx = 16;
constexpr int blocky = 8;

template <typename T>
void DotKillT::calcDiffMetric(PVideoFrame& f1, PVideoFrame& f2, int64_t* bdiffs, int nxblocks, int nyblocks, int field)
{
    const int plane[3] = { PLANAR_Y, PLANAR_U, PLANAR_V };
    const int planecount = std::min(vi.NumComponents(), 3);

    for (int i = 0; i < planecount; ++i)
    {
        const int stride = f1->GetPitch(plane[i]) / sizeof(T);
        const int width = f1->GetRowSize(plane[i]) / sizeof(T);
        const int height = f1->GetHeight(plane[i]);
        const T* f1p = reinterpret_cast<const T*>(f1->GetReadPtr(plane[i]));
        const T* f2p = reinterpret_cast<const T*>(f2->GetReadPtr(plane[i]));

        if (field)
        {
            f1p += stride;
            f2p += stride;
        }

        const int hblockx = (!i) ? blockx / 2 : (blockx / 2) / (1 << vi.GetPlaneWidthSubsampling(plane[i]));
        const int hblocky = (!i) ? blocky / 2 : (blocky / 2) / (1 << vi.GetPlaneHeightSubsampling(plane[i]));

        for (int y = 0; y < height / 2; ++y)
        {
            int ydest = y / hblocky;
            int xdest = 0;

            for (int x = 0; x < width; x += hblockx)
            {
                int acc = 0;
                int m = std::min(width, x + hblockx);
                for (int xl = x; xl < m; ++xl)
                {
                    int tmp = f1p[xl] - f2p[xl];
                    acc += tmp * tmp;
                }

                bdiffs[ydest * nxblocks + xdest] += acc;
                ++xdest;
            }

            f1p += stride * 2LL;
            f2p += stride * 2LL;
        }
    }
}

static int64_t getMaxDiff(int i, int j, const int64_t* bdiffs1, int nxblocks, int nyblocks)
{
    uint64_t tmp1 = bdiffs1[i * nxblocks + j] + bdiffs1[i * nxblocks + j + 1] + bdiffs1[(i + 1) * nxblocks + j] + bdiffs1[(i + 1) * nxblocks + j + 1];
    uint64_t tmp2 = bdiffs1[i * nxblocks + j] + bdiffs1[i * nxblocks + j - 1] + bdiffs1[(i + 1) * nxblocks + j] + bdiffs1[(i + 1) * nxblocks + j - 1];
    uint64_t tmp3 = bdiffs1[i * nxblocks + j] + bdiffs1[i * nxblocks - j + 1] + bdiffs1[(i - 1) * nxblocks + j] + bdiffs1[(i - 1) * nxblocks + j + 1];
    uint64_t tmp4 = bdiffs1[i * nxblocks + j] + bdiffs1[i * nxblocks - j - 1] + bdiffs1[(i - 1) * nxblocks + j] + bdiffs1[(i - 1) * nxblocks + j - 1];

    return std::max({ tmp1, tmp2, tmp3, tmp4 });
}

template <typename T>
static void diffMetricToMask(uint16_t* mask, const int64_t* bdiffs1, const int64_t* bdiffs2, int nxblocks, int nyblocks, int dupthresh, int tratio)
{
    int totdiff1 = 0;
    int totdiff2 = 0;

    for (int i = 1; i < nyblocks - 1; ++i)
    {
        for (int j = 1; j < nxblocks - 1; ++j)
        {
            uint64_t diff1 = getMaxDiff(i, j, bdiffs1, nxblocks, nyblocks);
            uint64_t diff2 = getMaxDiff(i, j, bdiffs2, nxblocks, nyblocks);

            if (diff1 >= dupthresh)
                ++totdiff1;
            if (diff2 >= dupthresh)
                ++totdiff2;
        }
    }

    // skip temporal processing if more than 1/tratio blocks have changed
    bool skip1 = (totdiff1 * tratio > (nxblocks - 2) * (nyblocks - 2));
    bool skip2 = (totdiff2 * tratio > (nxblocks - 2) * (nyblocks - 2));

    for (int i = 1; i < nyblocks - 1; ++i)
    {
        for (int j = 1; j < nxblocks - 1; ++j)
        {
            uint64_t diff1 = getMaxDiff(i, j, bdiffs1, nxblocks, nyblocks);
            uint64_t diff2 = getMaxDiff(i, j, bdiffs2, nxblocks, nyblocks);

            if (!skip1 && diff1 <= diff2 && diff1 < dupthresh)
                mask[nxblocks * i + j] = 1;
            else if (!skip2 && diff2 < diff1 && diff2 < dupthresh)
                mask[nxblocks * i + j] = 2;
            else
                mask[nxblocks * i + j] = 0;
        }

        // extend mask left and right
        mask[nxblocks * i] = mask[nxblocks * i + 1];
        mask[nxblocks * i + (nxblocks - 1)] = mask[nxblocks * i + (nxblocks - 2)];
    }

    // extend mask to top and bottom
    memcpy(mask, mask + nxblocks, nxblocks * sizeof(T));
    memcpy(mask + nxblocks * (nyblocks - 1LL), mask + nxblocks * (nyblocks - 2LL), nxblocks * sizeof(T));
}

template <typename T>
void DotKillT::applyTemporalMask(PVideoFrame& dst, PVideoFrame& f0, PVideoFrame& f1, PVideoFrame& f2, uint16_t* mask, int nxblocks, int nyblocks, int field, bool show)
{
    const int plane[3] = { PLANAR_Y, PLANAR_U, PLANAR_V };
    const int planecount = std::min(vi.NumComponents(), 3);

    for (int i = 0; i < planecount; ++i)
    {
        ptrdiff_t stride = f1->GetPitch(plane[i]) / sizeof(T);
        const int dst_stride = dst->GetPitch(plane[i]) / sizeof(T);
        const int width = f1->GetRowSize(plane[i]) / sizeof(T);
        const int height = f1->GetHeight(plane[i]);
        const T* f0p = reinterpret_cast<const T*>(f0->GetReadPtr(plane[i]));
        const T* f1p = reinterpret_cast<const T*>(f1->GetReadPtr(plane[i]));
        const T* f2p = reinterpret_cast<const T*>(f2->GetReadPtr(plane[i]));
        T* dstp = reinterpret_cast<T*>(dst->GetWritePtr(plane[i]));

        if (field)
        {
            f0p += stride;
            f1p += stride;
            f2p += stride;
            dstp += dst_stride;
        }

        const int hblockx = (!i) ? blockx / 2 : (blockx / 2) / (1 << vi.GetPlaneWidthSubsampling(plane[i]));
        const int hblocky = (!i) ? blocky / 2 : (blocky / 2) / (1 << vi.GetPlaneHeightSubsampling(plane[i]));

        for (int y = 0; y < height / 2; ++y)
        {
            int ydest = y / hblocky;
            int xdest = 0;

            for (int x = 0; x < width; x += hblockx)
            {
                int m = std::min(width, x + hblockx);

                for (int xl = x; xl < m; ++xl)
                {
                    if (mask[ydest * nxblocks + xdest] == 1)
                        dstp[xl] = ((f1p[xl] + f0p[xl] + 1) / 2);
                    else if (mask[ydest * nxblocks + xdest] == 2)
                        dstp[xl] = ((f2p[xl] + f0p[xl] + 1) / 2);
                    else
                        dstp[xl] = dstp[xl];
                }
                ++xdest;
            }

            f0p += stride * 2;
            f1p += stride * 2;
            f2p += stride * 2;
            dstp += dst_stride * 2;
        }
    }

    const int peak = (1 << vi.BitsPerComponent()) - 1;

    // Horrible square drawing code
    if (show)
    {
        ptrdiff_t stride = dst->GetPitch() / sizeof(T);
        T* dstp = reinterpret_cast<T*>(dst->GetWritePtr());

        const int width = dst->GetRowSize() / sizeof(T);
        const int height = dst->GetHeight();
        int hblockx = blockx / 2;

        for (int y = 0; y < height; ++y)
        {
            int ydest = y / blocky;
            int xdest = 0;

            for (int x = 0; x < width; x += hblockx)
            {
                int m = std::min(width, x + hblockx);

                if (y % blocky == 0)
                {
                    for (int xl = x; xl < m; ++xl)
                    {
                        if (mask[ydest * nxblocks + xdest] == 1)
                            dstp[xl] = 0;
                        else if (mask[ydest * nxblocks + xdest] == 2)
                            dstp[xl] = peak;
                    }
                }
                else if (y % blocky == blocky - 1)
                {
                    for (int xl = x; xl < m; ++xl)
                    {
                        if (mask[ydest * nxblocks + xdest] == 1)
                            dstp[xl] = 0;
                        else if (mask[ydest * nxblocks + xdest] == 2)
                            dstp[xl] = peak;
                    }
                }

                if (mask[ydest * nxblocks + xdest] == 1)
                {
                    dstp[x] = 0;
                    dstp[m - 1] = 0;
                }
                else if (mask[ydest * nxblocks + xdest] == 2)
                {
                    dstp[x] = 255;
                    dstp[m - 1] = peak;
                }

                ++xdest;
            }

            dstp += stride;
        }
    }
}

DotKillT::DotKillT(PClip _child, int _order, int _offset, int _dupthresh, int _tratio, bool _show, IScriptEnvironment* env)
    : GenericVideoFilter(_child), order(_order), offset(_offset), dupthresh(_dupthresh), tratio(_tratio), show(_show)
{
    if (vi.BitsPerComponent() == 32 || vi.IsRGB() || !vi.IsPlanar())
        env->ThrowError("DotKillT: only YUV 8..16-bit planar format supported.");
    if (order < 0 || order > 1)
        env->ThrowError("DotKillT: order must be either 0 or 1.");
    if (offset < 0 || offset > 4)
        env->ThrowError("DotKillT: offset must be between 0..4.");
    if (dupthresh < 0 || dupthresh > 255)
        env->ThrowError("DotKillT: dupthresh must be between 0..255.");

    has_at_least_v8 = true;
    try { env->CheckVersion(8); }
    catch (const AvisynthError&) { has_at_least_v8 = false; };

    if (vi.ComponentSize() == 1)
    {
        dupthresh *= dupthresh;

        applyfb = applyFieldBlend<uint8_t>;
        applydi = applyDotcrawInverse<uint8_t>;
        calcdif = &DotKillT::calcDiffMetric<uint8_t>;
        diffmetric = diffMetricToMask<uint8_t>;
        applytm = &DotKillT::applyTemporalMask<uint8_t>;
    }
    else
    {   
        dupthresh *= ((1 << vi.BitsPerComponent()) - 1) / 255;
        dupthresh *= dupthresh;

        applyfb = applyFieldBlend<uint16_t>;
        applydi = applyDotcrawInverse<uint16_t>;
        calcdif = &DotKillT::calcDiffMetric<uint16_t>;
        diffmetric = diffMetricToMask<uint16_t>;
        applytm = &DotKillT::applyTemporalMask<uint16_t>;
    }
}

PVideoFrame __stdcall DotKillT::GetFrame(int n, IScriptEnvironment* env)
{
    PVideoFrame srcpp = child->GetFrame(std::max(n - 2, 0), env);
    PVideoFrame srcp = child->GetFrame(std::max(n - 1, 0), env);
    PVideoFrame srcc = child->GetFrame(n, env);
    PVideoFrame srcn = child->GetFrame(std::min(n + 1, vi.num_frames - 1), env);
    PVideoFrame srcnn = child->GetFrame(std::min(n + 2, vi.num_frames - 1), env);

    // first two fields are duplicates, meaning that offset 0, 1 and 5, 6 are used in the various calculations
       // fields 2, 3 and 4 needs to have determined how many blocks are consecutively static
       // note that comparisons are only run on a single since we can in most comparisons can eliminate
       // the dotcrawl from the equation

       // 0-1 2 3 4 5-6

       // the complementary field can likewise be used for movement detection

       /*
           FIELD OFFSETS
            +  -  +  -  +  -
           -1  0  1  2  3  4
           A1 B1 B1 C1 D1 E1
           A2 B2 C2 D2 D2 E2
       */

    const int nxblocks = (vi.width + blockx / 2 - 1) / (blockx / 2);
    const int nyblocks = (vi.height + blocky / 2 - 1) / (blocky / 2);
    const int planecount = std::min(vi.NumComponents(), 3);

    PVideoFrame dst = srcc;
    env->MakeWritable(&dst);

    if (has_at_least_v8)
        env->propSetInt(env->getFramePropsRW(srcc), "DotKillTOffset", (n + offset) % 5, 0);

    if ((n + offset) % 5 == 0) {
        // 1
        applydi(srcc, srcn, dst, order, vi.BitsPerComponent(), planecount);

        // 2
        std::vector<int64_t> maskprev1(static_cast<int64_t>(nxblocks) * nyblocks * vi.ComponentSize());
        std::vector<int64_t> masknext1(static_cast<int64_t>(nxblocks) * nyblocks * vi.ComponentSize());
        std::vector<uint16_t> mask(static_cast<int64_t>(nxblocks) * nyblocks * vi.ComponentSize());
        (this->*calcdif)(srcp, srcn, maskprev1.data(), nxblocks, nyblocks, order);
        (this->*calcdif)(srcc, srcnn, masknext1.data(), nxblocks, nyblocks, order);

        diffmetric(mask.data(), maskprev1.data(), masknext1.data(), nxblocks, nyblocks, dupthresh, tratio);

        (this->*applytm)(dst, srcc, srcp, srcn, mask.data(), nxblocks, nyblocks, !order, show);
    }
    else if ((n + offset) % 5 == 1) {
        // 1
        applyfb(srcc, srcp, dst, order, planecount);

        // 2
        std::vector<int64_t> maskprev1(static_cast<int64_t>(nxblocks) * nyblocks * vi.ComponentSize());
        std::vector<int64_t> masknext1(static_cast<int64_t>(nxblocks) * nyblocks * vi.ComponentSize());
        std::vector<uint16_t> mask(static_cast<int64_t>(nxblocks) * nyblocks * vi.ComponentSize());
        (this->*calcdif)(srcp, srcn, maskprev1.data(), nxblocks, nyblocks, order);
        (this->*calcdif)(srcc, srcnn, masknext1.data(), nxblocks, nyblocks, !order);

        diffmetric(mask.data(), maskprev1.data(), masknext1.data(), nxblocks, nyblocks, dupthresh, tratio);

        (this->*applytm)(dst, srcc, srcp, srcn, mask.data(), nxblocks, nyblocks, !order, show);
    }
    else if ((n + offset) % 5 == 2) {
        // 1
        std::vector<int64_t> maskprev1(static_cast<int64_t>(nxblocks) * nyblocks * vi.ComponentSize());
        std::vector<int64_t> masknext1(static_cast<int64_t>(nxblocks) * nyblocks * vi.ComponentSize());
        std::vector<uint16_t> mask(static_cast<int64_t>(nxblocks) * nyblocks * vi.ComponentSize());
        (this->*calcdif)(srcc, srcpp, maskprev1.data(), nxblocks, nyblocks, order);
        (this->*calcdif)(srcp, srcn, masknext1.data(), nxblocks, nyblocks, !order);

        diffmetric(mask.data(), maskprev1.data(), masknext1.data(), nxblocks, nyblocks, dupthresh, tratio);

        (this->*applytm)(dst, srcc, srcp, srcn, mask.data(), nxblocks, nyblocks, order, show);

        // 2
        applyfb(srcc, srcn, dst, !order, planecount);
    }
    else if ((n + offset) % 5 == 3) {
        // 2
        applydi(srcc, srcp, dst, !order, vi.BitsPerComponent(), planecount);

        // 1
        std::vector<int64_t> maskprev1(static_cast<int64_t>(nxblocks) * nyblocks * vi.ComponentSize());
        std::vector<int64_t> masknext1(static_cast<int64_t>(nxblocks) * nyblocks * vi.ComponentSize());
        std::vector<uint16_t> mask(static_cast<int64_t>(nxblocks) * nyblocks * vi.ComponentSize());
        (this->*calcdif)(srcc, srcpp, maskprev1.data(), nxblocks, nyblocks, !order);
        (this->*calcdif)(srcp, srcn, masknext1.data(), nxblocks, nyblocks, !order);

        diffmetric(mask.data(), maskprev1.data(), masknext1.data(), nxblocks, nyblocks, dupthresh, tratio);

        (this->*applytm)(dst, srcc, srcp, srcn, mask.data(), nxblocks, nyblocks, order, show);
    }
    else if ((n + offset) % 5 == 4) {
        // 1
        {
            std::vector<int64_t> maskprev1(static_cast<int64_t>(nxblocks) * nyblocks * vi.ComponentSize());
            std::vector<int64_t> masknext1(static_cast<int64_t>(nxblocks) * nyblocks * vi.ComponentSize());
            std::vector<uint16_t> mask(static_cast<int64_t>(nxblocks) * nyblocks * vi.ComponentSize());
            (this->*calcdif)(srcc, srcpp, maskprev1.data(), nxblocks, nyblocks, !order);
            (this->*calcdif)(srcc, srcnn, masknext1.data(), nxblocks, nyblocks, order);

            diffmetric(mask.data(), maskprev1.data(), masknext1.data(), nxblocks, nyblocks, dupthresh, tratio);

            (this->*applytm)(dst, srcc, srcp, srcn, mask.data(), nxblocks, nyblocks, order, show);
        }

        // 2
        {
            std::vector<int64_t> maskprev1(static_cast<int64_t>(nxblocks) * nyblocks * vi.ComponentSize());
            std::vector<int64_t> masknext1(static_cast<int64_t>(nxblocks) * nyblocks * vi.ComponentSize());
            std::vector<uint16_t> mask(static_cast<int64_t>(nxblocks) * nyblocks * vi.ComponentSize());
            (this->*calcdif)(srcpp, srcc, maskprev1.data(), nxblocks, nyblocks, !order);
            (this->*calcdif)(srcc, srcnn, masknext1.data(), nxblocks, nyblocks, order);

            diffmetric(mask.data(), maskprev1.data(), masknext1.data(), nxblocks, nyblocks, dupthresh, tratio);

            (this->*applytm)(dst, srcc, srcp, srcn, mask.data(), nxblocks, nyblocks, !order, show);
        }
    }

    return dst;
}
