/**
 * @file   DumpUtils_test.cc
 * @brief  Unit test for utilities in DumpUtils.h
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   February 28, 2017
 *
 *
 */

// Boost libraries
#define BOOST_TEST_MODULE DumpUtils_test
#include <boost/test/unit_test.hpp>

// LArSoft libraries
#include "larcorealg/CoreUtils/DumpUtils.h"

// ROOT libraries
#include <TVector3.h>

// C/C++ standard libraries
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <array>
#include <sstream>


//------------------------------------------------------------------------------
void ArrayTest_STLarray() {
  std::ostringstream sstr;
  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  std::array<float, 5> const array = {{ 1., 2., 3., 4., 6. }};
  sstr << lar::dump::array<5>(array);
  BOOST_CHECK_EQUAL(sstr.str(), "{ 1; 2; 3; 4; 6 }");
  BOOST_CHECK_EQUAL
    (std::string(lar::dump::array<5>(array)), "{ 1; 2; 3; 4; 6 }");
}

void ArrayTest_STLvector() {
  std::ostringstream sstr;
  std::vector<float> const array = { 1., 2., 3., 4., 6. };
  sstr << lar::dump::array<5>(array);
  BOOST_CHECK_EQUAL(sstr.str(), "{ 1; 2; 3; 4; 6 }");
  BOOST_CHECK_EQUAL
    (std::string(lar::dump::array<5>(array)), "{ 1; 2; 3; 4; 6 }");
}

void ArrayTest_carray() {
  std::ostringstream sstr;
  float const array[5] = { 1., 2., 3., 4., 6. };
  sstr << lar::dump::array<5>(array);
  BOOST_CHECK_EQUAL(sstr.str(), "{ 1; 2; 3; 4; 6 }");
  BOOST_CHECK_EQUAL
    (std::string(lar::dump::array<5>(array)), "{ 1; 2; 3; 4; 6 }");
}

void ArrayTest_cptr() {
  std::ostringstream sstr;
  float const array[5] = { 1., 2., 3., 4., 6. };
  float const* ptr = array;
  sstr << lar::dump::array<5>(ptr);
  BOOST_CHECK_EQUAL(sstr.str(), "{ 1; 2; 3; 4; 6 }");
  BOOST_CHECK_EQUAL(std::string(lar::dump::array<5>(ptr)), "{ 1; 2; 3; 4; 6 }");
}

void ArrayTest_array() {
  std::ostringstream sstr;
  double array[5] = { 1., 2., 3., 4., 6. };
  sstr << lar::dump::array<5>(array);
  BOOST_CHECK_EQUAL(sstr.str(), "{ 1; 2; 3; 4; 6 }");
  BOOST_CHECK_EQUAL
    (std::string(lar::dump::array<5>(array)), "{ 1; 2; 3; 4; 6 }");
}

void ArrayTest_ptr() {
  std::ostringstream sstr;
  double array[5] = { 1., 2., 3., 4., 6. };
  double* ptr = array;
  sstr << lar::dump::array<5>(ptr);
  BOOST_CHECK_EQUAL(sstr.str(), "{ 1; 2; 3; 4; 6 }");
  BOOST_CHECK_EQUAL(std::string(lar::dump::array<5>(ptr)), "{ 1; 2; 3; 4; 6 }");
}

void ArrayTest() {

  ArrayTest_STLarray();
  ArrayTest_STLvector();
  ArrayTest_array();
  ArrayTest_carray();
  ArrayTest_ptr();
  ArrayTest_cptr();

} // ArrayTest()


//------------------------------------------------------------------------------
void ArrayDocumentationTest() {

  std::ostringstream sstr;

  /* This is the code as in the documentation:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * double const data[5] = { 1., 2., 3., 4., 6. };
   * std::cout << "Data: " << lar::dump::array<5>(data) << std::endl;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * which is promised to return:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * Data: { 1; 2; 3; 4; 6 }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */

  double const data[5] = { 1., 2., 3., 4., 6. };
  std::cout << "Data: " << lar::dump::array<5>(data) << std::endl;

  sstr << lar::dump::array<5>(data);
  BOOST_CHECK_EQUAL(sstr.str(), "{ 1; 2; 3; 4; 6 }");

} // ArrayDocumentationTest()


//------------------------------------------------------------------------------
void VectorTest_STLvector() {
  std::ostringstream sstr;
  std::vector<float> const v = { 1., 2., 3., 4., 6. };
  sstr << lar::dump::vector(v);
  BOOST_CHECK_EQUAL(sstr.str(), "{ 1; 2; 3; 4; 6 }");
  BOOST_CHECK_EQUAL
    (std::string(lar::dump::vector(v)), "{ 1; 2; 3; 4; 6 }");
} // VectorTest_STLvector()

void VectorTest_STLlist() {
  std::ostringstream sstr;
  std::list<float> const v = { 1., 2., 3., 4., 6. };
  sstr << lar::dump::vector(v);
  BOOST_CHECK_EQUAL(sstr.str(), "{ 1; 2; 3; 4; 6 }");
  BOOST_CHECK_EQUAL
    (std::string(lar::dump::vector(v)), "{ 1; 2; 3; 4; 6 }");
} // VectorTest_STLlist()

void VectorTest() {

  VectorTest_STLvector();
  VectorTest_STLlist();

} // VectorTest()


//------------------------------------------------------------------------------
void VectorDocumentationTest() {

  std::ostringstream sstr;

  /* This is the code as in the documentation:
   *
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<double> data = { 1., 2., 3., 4., 6. };
   * std::cout << "Data: " << lar::dump::vector(data) << std::endl;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * which is promised to return:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * Data: { 1; 2; 3; 4; 6 }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */

  std::vector<double> data = { 1., 2., 3., 4., 6. };
  std::cout << "Data: " << lar::dump::vector(data) << std::endl;

  sstr << lar::dump::array<5>(data);
  BOOST_CHECK_EQUAL(sstr.str(), "{ 1; 2; 3; 4; 6 }");

} // VectorDocumentationTest()


//------------------------------------------------------------------------------
void Vector3Dtest_TVector3() {
  std::ostringstream sstr;
  TVector3 pos( 1., 2., 4. );
  sstr << lar::dump::vector3D(pos);
  BOOST_CHECK_EQUAL(sstr.str(), "{ 1; 2; 4 }");
  BOOST_CHECK_EQUAL(std::string(lar::dump::vector3D(pos)), "{ 1; 2; 4 }");
}

void Vector3Dtest_carray() {
  std::ostringstream sstr;
  float const array[3] = { 1., 2., 4. };
  sstr << lar::dump::vector3D(array);
  BOOST_CHECK_EQUAL(sstr.str(), "{ 1; 2; 4 }");
  BOOST_CHECK_EQUAL(std::string(lar::dump::vector3D(array)), "{ 1; 2; 4 }");
}

void Vector3Dtest_cptr() {
  std::ostringstream sstr;
  float const array[3] = { 1., 2., 4. };
  float const* ptr = array;
  sstr << lar::dump::vector3D(ptr);
  BOOST_CHECK_EQUAL(sstr.str(), "{ 1; 2; 4 }");
  BOOST_CHECK_EQUAL(std::string(lar::dump::vector3D(ptr)), "{ 1; 2; 4 }");
}

void Vector3Dtest_array() {
  std::ostringstream sstr;
  double array[3] = { 1., 2., 4. };
  sstr << lar::dump::vector3D(array);
  BOOST_CHECK_EQUAL(sstr.str(), "{ 1; 2; 4 }");
  BOOST_CHECK_EQUAL(std::string(lar::dump::vector3D(array)), "{ 1; 2; 4 }");
}

void Vector3Dtest_ptr() {
  std::ostringstream sstr;
  double array[3] = { 1., 2., 4. };
  double* ptr = array;
  sstr << lar::dump::vector3D(ptr);
  BOOST_CHECK_EQUAL(sstr.str(), "{ 1; 2; 4 }");
  BOOST_CHECK_EQUAL(std::string(lar::dump::vector3D(ptr)), "{ 1; 2; 4 }");
}


void Vector3Dtest() {

  Vector3Dtest_TVector3();
  Vector3Dtest_array();
  Vector3Dtest_carray();
  Vector3Dtest_ptr();
  Vector3Dtest_cptr();

} // Vector3Dtest()


//------------------------------------------------------------------------------
void Vector3DstreamOutputDocumentationTest() {

  std::ostringstream sstr;

  /* This is the code as in the documentation:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * TVector3 pos;
   * std::cout << "Position: " << lar::dump::vector3D(pos) << std::endl;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * which is promised to return:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * Position: { 0, 0, 0 }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */

  TVector3 pos;
  std::cout << "Position: " << lar::dump::vector3D(pos) << std::endl;

  sstr << lar::dump::vector3D(pos);
  BOOST_CHECK_EQUAL(sstr.str(), "{ 0; 0; 0 }");

} // Vector3DstreamOutputDocumentationTest()


void Vector3DstringConcatDocumentationTest() {

  std::ostringstream sstr;

  /* This is the code as in the documentation:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * TVector3 pos;
   * std::runtime_error e("Position: " + lar::dump::vector3D(pos));
   * std::cout << e.what() << std::endl;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * which is promised to return:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * Position: { 0, 0, 0 }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */

  TVector3 pos;
  std::runtime_error e("Position: " + lar::dump::vector3D(pos));
  std::cout << e.what() << std::endl;

  BOOST_CHECK_EQUAL(std::string(e.what()), "Position: { 0; 0; 0 }");

} // Vector3DstringConcatDocumentationTest()


void Vector3DstringAppendDocumentationTest() {

  std::ostringstream sstr;

  /* This is the code that should have been in the documentation:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * TVector3 pos;
   * std::string s = "Position: ";
   * s += lar::dump::vector3D(pos);
   * std::cout << s << std::endl;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * which is promised to return:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * Position: { 0, 0, 0 }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */

  TVector3 pos;
  std::string s = "Position: ";
  s += lar::dump::vector3D(pos);
  std::cout << s << std::endl;

  BOOST_CHECK_EQUAL(s, "Position: { 0; 0; 0 }");

} // Vector3DstringAppendDocumentationTest()


void Vector3DdocumentationTest() {
  Vector3DstreamOutputDocumentationTest();
  Vector3DstringConcatDocumentationTest();
  Vector3DstringAppendDocumentationTest();
} // Vector3DdocumentationTest()


//------------------------------------------------------------------------------

#include <array>

namespace lar {
  namespace dump {

    // specialization for 3D std::array<>
    template <typename T>
    struct VectorDumper<std::array<T, 3U>> {

      using Vector_t = std::array<T, 3U>;

      Vector_t const& a;

      VectorDumper(Vector_t const& a): a(a) {}

      template <typename Stream>
      void operator() (Stream&& out) const
        { out << "{ " << a[0] << "; " << a[1] << "; " << a[2] << " }"; }

    }; // struct VectorDumper<array>


  } // namespace dump
} // namespace lar


void Vector3DspecializationTest() {

  std::ostringstream sstr;

  std::array<float, 3U> v;
  v.fill(0.0);
  sstr << lar::dump::vector3D(v);
  BOOST_CHECK_EQUAL(sstr.str(), "{ 0; 0; 0 }");

} // Vector3DspecializationTest()



//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(DumpArrays_testcase) {

  ArrayTest();

  ArrayDocumentationTest();

} // BOOST_AUTO_TEST_CASE(DumpArrays_testcase)


BOOST_AUTO_TEST_CASE(DumpVectors_testcase) {

  VectorTest();

  VectorDocumentationTest();

} // BOOST_AUTO_TEST_CASE(DumpVectors_testcase)


BOOST_AUTO_TEST_CASE(Dump3Dvectors_testcase) {

  Vector3Dtest();

  Vector3DdocumentationTest();

  Vector3DspecializationTest();

} // BOOST_AUTO_TEST_CASE(Dump3Dvectors_testcase)

