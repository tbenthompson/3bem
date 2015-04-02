#ifndef __UNITTEST_SHARED_H
#define __UNITTEST_SHARED_H

#include "UnitTest++.h"
#include "TestReporterStdout.h"
#include <cstring>

int RunOneTest(std::string name) {
    UnitTest::TestReporterStdout reporter;
    UnitTest::TestRunner runner(reporter);
    return runner.RunTestsIf(UnitTest::Test::GetTestList(), nullptr, 
        [&](const UnitTest::Test* const test) {
            return (0 == strcmp(test->m_details.testName, name.c_str()));
        }, 0);
}
#endif
