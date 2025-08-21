#include <cstdio>
#include <cstdlib>
#include <string>
#include "uSignal/version.hpp"

using namespace USignal;

int Version::getMajor() noexcept
{
    return USIGNAL_MAJOR;
}

int Version::getMinor() noexcept
{
    return USIGNAL_MINOR;
}

int Version::getPatch() noexcept
{
    return USIGNAL_PATCH;
}

bool Version::isAtLeast(const int major, const int minor,
                        const int patch) noexcept
{
    if (USIGNAL_MAJOR < major){return false;}
    if (USIGNAL_MAJOR > major){return true;}
    if (USIGNAL_MINOR < minor){return false;}
    if (USIGNAL_MINOR > minor){return true;}
    if (USIGNAL_PATCH < patch){return false;}
    return true;
}

std::string Version::getVersion() noexcept
{
    std::string version(USIGNAL_VERSION);
    return version;
}
