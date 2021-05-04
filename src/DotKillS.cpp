#include <cstring>

#include "DotKill.h"

template <typename T>
static void convHoriz(const uint8_t* src_, int srcStride, int* dst, int width, int height)
{
    const T* src = reinterpret_cast<const T*>(src_);

    while (--height)
    {
        dst[0] = 0;
        dst[1] = 0;

        for (int x = 2; x < width - 3; ++x)
        {
            int temp = -(src[x - 2] + src[x - 1]) + 2 * (src[x] + src[x + 1]) - (src[x + 2] + src[x + 3]);
            temp += -(src[x - 2 + srcStride] + src[x - 1 + srcStride]) + 2 * (src[x + srcStride] + src[x + 1 + srcStride]) - (src[x + 2 + srcStride] + src[x + 3 + srcStride]);
            dst[x] = temp;
        }

        dst[width - 3] = 0;
        dst[width - 2] = 0;
        dst[width - 1] = 0;

        src += srcStride;
        dst += width;
    }

    memset(dst, 0, sizeof(int) * width);
}

template <typename T>
static void convVert(const uint8_t* src_, int srcStride, int* dst, int width, int height)
{
    height -= 5;
    const T* src = reinterpret_cast<const T*>(src_);

    memset(dst, 0, sizeof(int) * width * 2);
    src += 2LL * srcStride;
    dst += 2LL * width;

    while (height--)
    {
        for (int x = 0; x < width - 1; ++x)
        {
            dst[x] = -(src[x - 2 * srcStride] + src[x - 2 * srcStride + 1] + (src[x - 1 * srcStride] + src[x - 1 * srcStride + 1]))
                + 2 * (src[x] + src[x + 1] + (src[x + 1 * srcStride] + src[x + 1 * srcStride + 1]))
                - (src[x + 2 * srcStride] + src[x + 2 * srcStride + 1] + (src[x + 3 * srcStride] + src[x + 3 * srcStride + 1]));
        }

        dst[width - 1] = 0;

        src += srcStride;
        dst += width;
    }

    memset(dst, 0, sizeof(int) * width * 3);
}

template <typename T>
static void applyMask(const int* maskPtr, uint8_t* dst_, int dstStride, int width, int height, uint8_t* ppMask_, int depth)
{
    T* dst = reinterpret_cast<T*>(dst_);
    T* ppMask = reinterpret_cast<T*>(ppMask_);

    maskPtr += width;
    ppMask += width;
    dst += dstStride;

    const int shift = depth - 8;
    const int lmax = 235 << shift;
    const int lmin = 16 << shift;
    const int peak = (1 << depth) - 1;
    int sortArray[8];

    for (int y = 1; y < height - 1; ++y)
    {
        for (int x = 1; x < width - 2; ++x)
        {
            sortArray[0] = maskPtr[x - width - 1];
            sortArray[1] = maskPtr[x - width];
            sortArray[2] = maskPtr[x - width + 1];
            sortArray[3] = maskPtr[x - 1];
            sortArray[4] = maskPtr[x + 1];
            sortArray[5] = maskPtr[x + width - 1];
            sortArray[6] = maskPtr[x + width];
            sortArray[7] = maskPtr[x + width + 1];

            std::sort(sortArray, sortArray + 8);

            int upper = sortArray[7];
            int lower = sortArray[0];

            int t = maskPtr[x] - std::clamp(maskPtr[x], lower, upper);

            if (t >= 0)
                t = (t + 4) / 8;
            else
                t = (t - 4) / 8;

            if (std::abs(t) > 1)
            {
                ppMask[x] = peak;

                dst[x] = static_cast<T>(std::clamp(dst[x] - t, lmin, lmax));
                dst[x + 1] = static_cast<T>(std::clamp(dst[x + 1] - t, lmin, lmax));
                dst[x + dstStride] = static_cast<T>(std::clamp(dst[x + dstStride] - t, lmin, lmax));
                dst[x + 1 + dstStride] = static_cast<T>(std::clamp(dst[x + 1 + dstStride] - t, lmin, lmax));
            }
        }

        maskPtr += width;
        ppMask += width;
        dst += dstStride;
    }
}

DotKillS::DotKillS(PClip _child, int _iterations, IScriptEnvironment* env)
    : GenericVideoFilter(_child), iterations(_iterations)
{
    if (vi.BitsPerComponent() == 32 || vi.IsRGB() || !vi.IsPlanar())
        env->ThrowError("DotKillS: only YUV 8..16-bit planar format supported.");
    if (iterations < 1)
        env->ThrowError("DotKillS: iterations must be equal to or greater than 1.");

    if (vi.ComponentSize() == 1)
    {
        ch = convHoriz<uint8_t>;
        cv = convVert<uint8_t>;
        apply = applyMask<uint8_t>;
    }
    else
    {
        ch = convHoriz<uint16_t>;
        cv = convVert<uint16_t>;
        apply = applyMask<uint16_t>;
    }
}

PVideoFrame __stdcall DotKillS::GetFrame(int n, IScriptEnvironment* env)
{
    PVideoFrame frame = child->GetFrame(n, env);
    env->MakeWritable(&frame);

    const int dst_stride = frame->GetPitch() / vi.ComponentSize();
    const int width = vi.width;
    const int height = vi.height;
    uint8_t* dstp = frame->GetWritePtr();

    int* tempMask = new int[static_cast<int64_t>(width) * height];
    uint8_t* ppMask = new uint8_t[static_cast<int64_t>(width) * height * vi.ComponentSize()]();

    for (int i = 0; i < iterations; ++i)
    {
        cv(dstp, dst_stride, tempMask, width, height);
        apply(tempMask, dstp, dst_stride, width, height, ppMask, vi.BitsPerComponent());

        ch(dstp, dst_stride, tempMask, width, height);
        apply(tempMask, dstp, dst_stride, width, height, ppMask, vi.BitsPerComponent());
    }

    delete[] tempMask;
    delete[] ppMask;

    return frame;
}
