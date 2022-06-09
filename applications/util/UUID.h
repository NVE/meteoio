// https://stackoverflow.com/a/58467162

#ifndef UUID_H
#define UUID_H

#include <random>

using namespace std;

class UUID
{
public:
    static string generate()
    {
        static random_device dev;
        static mt19937 rng(dev());

        uniform_int_distribution<int> dist(0, 15);

        const char *v = "0123456789abcdef";
        const bool dash[] = {0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0};

        string uuid;
        for (int i = 0; i < 16; i++)
        {
            if (dash[i])
                uuid += "-";
            uuid += v[dist(rng)];
            uuid += v[dist(rng)];
        }
        return uuid;
    }

private:
    UUID() {}
};

#endif // UUID_H