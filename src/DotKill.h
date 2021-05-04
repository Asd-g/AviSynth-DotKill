#pragma once

#include <algorithm>

#include "avisynth.h"

class DotKillS : public GenericVideoFilter
{
    int iterations;

    void (*ch)(const uint8_t* src, int srcStride, int* dst, int width, int height);
    void (*cv)(const uint8_t* src, int srcStride, int* dst, int width, int height);
    void (*apply)(const int* maskPtr, uint8_t* dst, int dstStride, int width, int height, uint8_t* ppMask, int depth);

public:
    DotKillS(PClip _child, int _iterations, IScriptEnvironment* env);

    int __stdcall SetCacheHints(int cachehints, int frame_range)
    {
        return cachehints == CACHE_GET_MTMODE ? MT_NICE_FILTER : 0;
    }

    PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);
};

class DotKillZ : public GenericVideoFilter
{
    int order;
    int offset;

    void (*applyfb)(PVideoFrame& srcc, PVideoFrame& srcn, PVideoFrame& dst, int order, int planecount);
    void (*applydi)(PVideoFrame& srcc, PVideoFrame& srcn, PVideoFrame& dst, int order, int depth, int planecount);

public:
    DotKillZ(PClip _child, int _order, int _offset, IScriptEnvironment* env);

    int __stdcall SetCacheHints(int cachehints, int frame_range)
    {
        return cachehints == CACHE_GET_MTMODE ? MT_NICE_FILTER : 0;
    }

    PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);
};

class DotKillT : public GenericVideoFilter
{
    int order;
    int offset;
    int dupthresh;
    int tratio;
    bool show;
    bool has_at_least_v8;

    template <typename T>
    void calcDiffMetric(PVideoFrame& f1, PVideoFrame& f2, int64_t* bdiffs, int nxblocks, int nyblocks, int field);
    template <typename T>
    void applyTemporalMask(PVideoFrame& dst, PVideoFrame& f0, PVideoFrame& f1, PVideoFrame& f2, uint16_t* mask, int nxblocks, int nyblocks, int field, bool show);

    void (*applyfb)(PVideoFrame& srcc, PVideoFrame& srcn, PVideoFrame& dst, int order, int planecount);
    void (*applydi)(PVideoFrame& srcc, PVideoFrame& srcn, PVideoFrame& dst, int order, int depth, int planecount);
    void (DotKillT::* calcdif)(PVideoFrame& f1, PVideoFrame& f2, int64_t* bdiffs, int nxblocks, int nyblocks, int field);
    void (*diffmetric)(uint16_t* mask_, const int64_t* bdiffs1, const int64_t* bdiffs2, int nxblocks, int nyblocks, int dupthresh, int tratio);
    void (DotKillT::* applytm)(PVideoFrame& dst, PVideoFrame& f0, PVideoFrame& f1, PVideoFrame& f2, uint16_t* mask, int nxblocks, int nyblocks, int field, bool show);

public:
    DotKillT(PClip _child, int _order, int _offset, int _dupthresh, int _tratio, bool _show, IScriptEnvironment* env);

    int __stdcall SetCacheHints(int cachehints, int frame_range)
    {
        return cachehints == CACHE_GET_MTMODE ? MT_NICE_FILTER : 0;
    }

    PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);
};

template <typename T>
void applyFieldBlend(PVideoFrame& srcc, PVideoFrame& srcn, PVideoFrame& dst, int order, int planecount);
template <typename T>
void applyDotcrawInverse(PVideoFrame& srcc, PVideoFrame& srcn, PVideoFrame& dst, int order, int depth, int planecount);
