// Force-included (/FI) into the MVB++ sources on MSVC only.
//
// <windows.h> (reached through OCCT) defines the 16-bit leftovers `near` and `far` as empty macros, so
// ordinary identifiers such as `TopoDS_Shape far = ...` or `sp.near` stop compiling. Include it here
// first, under its own include guard, then drop the lowercase macros. NEAR/FAR stay defined (empty) so
// the Windows headers included later still parse.
//
// It is included lean: NOGDI keeps out wingdi.h, whose ERROR macro breaks MKF's LogLevel::ERROR (and
// MVB++ draws nothing through GDI); WIN32_LEAN_AND_MEAN drops the rarely used APIs; NOMINMAX keeps
// std::min/std::max usable.
#pragma once
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOGDI
#define NOGDI
#endif
#include <windows.h>
#undef near
#undef far
#undef NEAR
#undef FAR
#define NEAR
#define FAR
