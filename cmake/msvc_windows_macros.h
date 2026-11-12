// Force-included (/FI) into the MVB++ sources on MSVC only.
//
// <windows.h> (reached through OCCT) defines the 16-bit leftovers `near` and `far` as empty macros, so
// ordinary identifiers such as `TopoDS_Shape far = ...` or `sp.near` stop compiling. Include it here
// first, under its own include guard, then drop the lowercase macros. NEAR/FAR stay defined (empty) so
// the Windows headers included later still parse. NOMINMAX keeps std::min/std::max usable.
#pragma once
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#undef near
#undef far
#undef NEAR
#undef FAR
#define NEAR
#define FAR
