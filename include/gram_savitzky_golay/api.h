#pragma once

// Handle portable symbol export.
// Defining manually which symbol should be exported is required
// under Windows whether MinGW or MSVC is used.
//
// The headers then have to be able to work in two different modes:
// - dllexport when one is building the library,
// - dllimport for clients using the library.
//
// On Linux, set the visibility accordingly. If C++ symbol visibility
// is handled by the compiler, see: http://gcc.gnu.org/wiki/Visibility
#if defined _WIN32 || defined __CYGWIN__
// On Microsoft Windows, use dllimport and dllexport to tag symbols.
#  define GRAM_SAVITZKY_GOLAY_DLLIMPORT __declspec(dllimport)
#  define GRAM_SAVITZKY_GOLAY_DLLEXPORT __declspec(dllexport)
#  define GRAM_SAVITZKY_GOLAY_DLLLOCAL
#else
// On Linux, for GCC >= 4, tag symbols using GCC extension.
#  if __GNUC__ >= 4
#    define GRAM_SAVITZKY_GOLAY_DLLIMPORT __attribute__((visibility("default")))
#    define GRAM_SAVITZKY_GOLAY_DLLEXPORT __attribute__((visibility("default")))
#    define GRAM_SAVITZKY_GOLAY_DLLLOCAL __attribute__((visibility("hidden")))
#  else
// Otherwise (GCC < 4 or another compiler is used), export everything.
#    define GRAM_SAVITZKY_GOLAY_DLLIMPORT
#    define GRAM_SAVITZKY_GOLAY_DLLEXPORT
#    define GRAM_SAVITZKY_GOLAY_DLLLOCAL
#  endif // __GNUC__ >= 4
#endif // defined _WIN32 || defined __CYGWIN__

#ifdef GRAM_SAVITZKY_GOLAY_STATIC
// If one is using the library statically, get rid of
// extra information.
#  define GRAM_SAVITZKY_GOLAY_DLLAPI
#  define GRAM_SAVITZKY_GOLAY_LOCAL
#else
// Depending on whether one is building or using the
// library define DLLAPI to import or export.
#  ifdef GRAM_SAVITZKY_GOLAY_EXPORTS
#    define GRAM_SAVITZKY_GOLAY_DLLAPI GRAM_SAVITZKY_GOLAY_DLLEXPORT
#  else
#    define GRAM_SAVITZKY_GOLAY_DLLAPI GRAM_SAVITZKY_GOLAY_DLLIMPORT
#  endif // GRAM_SAVITZKY_GOLAY_EXPORTS
#  define GRAM_SAVITZKY_GOLAY_LOCAL GRAM_SAVITZKY_GOLAY_DLLLOCAL
#endif // GRAM_SAVITZKY_GOLAY_STATIC
