#include "gregfile.hpp"
#include <iostream>
#include <iomanip>
#include <cstring>

int main() {
	gtd::file f{"out.txt"};
	std::cout << std::boolalpha << f.write("hi", 2) << std::endl;;
	f.close();
    gtd::file r{"out.txt", O_RDONLY};
    char buff[8192]{};
    std::cout << std::boolalpha << r.read(buff) << std::endl;;
    perror("Error: ");
    r.close();
    std::cout << buff << std::endl;
	return 0;
}

