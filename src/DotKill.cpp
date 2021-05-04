#include "DotKill.h"

AVSValue __cdecl Create_DotKillS(AVSValue args, void* user_data, IScriptEnvironment* env)
{
    return new DotKillS(
        args[0].AsClip(),
        args[1].AsInt(1),
        env);
}

AVSValue __cdecl Create_DotKillZ(AVSValue args, void* user_data, IScriptEnvironment* env)
{
    return new DotKillZ(
        args[0].AsClip(),
        args[1].AsInt(0),
        args[2].AsInt(0),
        env);
}

AVSValue __cdecl Create_DotKillT(AVSValue args, void* user_data, IScriptEnvironment* env)
{
    return new DotKillT(
        args[0].AsClip(),
        args[1].AsInt(0),
        args[2].AsInt(0),
        args[3].AsInt(64),
        args[4].AsInt(3),
        args[5].AsBool(false),
        env);
}

const AVS_Linkage* AVS_linkage;

extern "C" __declspec(dllexport)
const char* __stdcall AvisynthPluginInit3(IScriptEnvironment * env, const AVS_Linkage* const vectors)
{
    AVS_linkage = vectors;

    env->AddFunction("DotKillS", "c[iterations]i", Create_DotKillS, 0);
    env->AddFunction("DotKillZ", "c[order]i[offset]i", Create_DotKillZ, 0);
    env->AddFunction("DotKillT", "c[order]i[offset]i[dupthresh]i[tratio]i[show]b", Create_DotKillT, 0);

    return "DotKill";
}
