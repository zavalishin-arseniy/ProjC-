# CMakeList.txt: ������ CMake ��� C++Test; �������� �������� ��� � �����������,
# ������� ����� ������ ��� ����������� �������.
#
#cmake_minimum_required (VERSION 3.8)

# �������� �������� � ����������� ���� ����� �������.
#add_executable (C++Test "C++Test.cpp" )

# TODO: �������� ����� � ������� �������, ���� ��� ����������.


cmake_minimum_required(VERSION 2.8)
project(ProjC-)

set(SOURCE_EXE main.cpp)

set(SOURCE_LIB calc.cpp)

add_library(calc STATIC ${SOURCE_LIB})

add_executable(main ${SOURCE_EXE})

add_executable(tests tests.cpp)

target_link_libraries(tests calc)

add_test(tests "./tests")

target_link_libraries(main calc)