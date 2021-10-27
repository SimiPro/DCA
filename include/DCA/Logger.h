#pragma once

#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include <string>
#include <vector>

namespace DCA {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#define ANSI_COLOR_RED "\x1b[31m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_YELLOW "\x1b[33m"
#define ANSI_COLOR_BLUE "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN "\x1b[36m"
#define ANSI_COLOR_DEFAULT "\x1b[0m"

#define GET_STRING_FROM_ARGUMENT_LIST(fmt, pBuffer)              \
    {                                                            \
        va_list args;                                            \
        pBuffer = nullptr;                                       \
        int length = 1024;                                       \
        int result = -1;                                         \
        while (result == -1) {                                   \
            delete[] pBuffer;                                    \
            pBuffer = new char[length + 1];                      \
            memset(pBuffer, 0, length + 1);                      \
            va_start(args, fmt);                                 \
            result = std::vsnprintf(pBuffer, length, fmt, args); \
            va_end(args);                                        \
            if (result >= length)                                \
                result = -1;                                     \
            length *= 2;                                         \
        }                                                        \
    }

#define RELEASE_STRING_FROM_ARGUMENT_LIST(pBuffer) \
    {                                              \
        delete[] pBuffer;                          \
        pBuffer = nullptr;                         \
    }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 * @brief %Logger class.
 * 
 * Writes to stout.
 */
class Logger {
public:
    /**
     * @brief Possible color values.
     */
    enum PRINT_COLOR { RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, DEFAULT };

    /**
     * @brief Print a message.
     * @param[in] color The color to use.
     * @param[in] fmt The thing to print, potentially containing format strings.
     */
    static void print(PRINT_COLOR color, const char* fmt, ...) {
        char* pBuffer = nullptr;
        GET_STRING_FROM_ARGUMENT_LIST(fmt, pBuffer);

        if (color == Logger::RED)
            printf("%s%s%s", ANSI_COLOR_RED, pBuffer, ANSI_COLOR_DEFAULT);
        else if (color == Logger::GREEN)
            printf("%s%s%s", ANSI_COLOR_GREEN, pBuffer, ANSI_COLOR_DEFAULT);
        else if (color == Logger::YELLOW)
            printf("%s%s%s", ANSI_COLOR_YELLOW, pBuffer, ANSI_COLOR_DEFAULT);
        else if (color == Logger::BLUE)
            printf("%s%s%s", ANSI_COLOR_BLUE, pBuffer, ANSI_COLOR_DEFAULT);
        else if (color == Logger::MAGENTA)
            printf("%s%s%s", ANSI_COLOR_MAGENTA, pBuffer, ANSI_COLOR_DEFAULT);
        else if (color == Logger::CYAN)
            printf("%s%s%s", ANSI_COLOR_CYAN, pBuffer, ANSI_COLOR_DEFAULT);
        else if (color == Logger::DEFAULT)
            printf("%s%s%s", ANSI_COLOR_DEFAULT, pBuffer, ANSI_COLOR_DEFAULT);

        RELEASE_STRING_FROM_ARGUMENT_LIST(pBuffer);
    }
};

}  // namespace DCA