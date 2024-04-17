#ifndef WORLD_HPP
#define WORLD_HPP

#include <string>

class World
{
public:
    void set(std::string msg) { this->message = msg; }
    std::string greet() { return message; }

private:
    std::string message;
};

#endif // WORLD_HPP