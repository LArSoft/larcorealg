#include "larcorealg/CoreUtils/DebugUtils.h"

// C/C++ standard libraries
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct AClass {
  void Do()
  {
    std::cout << present() << ": Do()" << std::endl;
    lar::debug::printBacktrace(std::cout);
  }

  std::string present() const
  {
    std::ostringstream sstr;
    sstr << lar::debug::demangle(this) << "[" << ((void*)this) << "]";
    return sstr.str();
  } // present()

}; // AClass

struct BClass {
  void Do()
  {
    std::cout << present() << ": Do()" << std::endl;
    lar::debug::printBacktrace(std::cout);
    a.Do();
  }
  std::string present() const
  {
    std::ostringstream sstr;
    sstr << lar::debug::demangle(this) << "[" << ((void*)this) << "]";
    return sstr.str();
  } // present()

  AClass a;
}; // BClass

struct CClass {
  CClass() { Do(); }

  void Do()
  {
    std::cout << present() << ": Do()" << std::endl;
    lar::debug::printBacktrace(std::cout);
    a.Do();
  }
  std::string present() const
  {
    std::ostringstream sstr;
    sstr << lar::debug::demangle(this) << "[" << ((void*)this) << "]";
    return sstr.str();
  } // present()

  AClass a;
}; // CClass

int main(int /* argc */, char** /* argv */)
{

  AClass a;
  a.Do();

  BClass b;
  b.Do();

  // create test classes in a deeper, non-local context
  std::vector<CClass> v [[maybe_unused]] (2);

  return 0;
} // main()
