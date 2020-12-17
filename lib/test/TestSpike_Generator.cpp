#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <cppunit/TestCase.h>
#include <cppunit/TestFixture.h>
#include <cppunit/ui/text/TextTestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>
#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/XmlOutputter.h>
#include <netinet/in.h>

#include "../include/spike_generator.hpp"

using namespace CppUnit;
using namespace std;

class TestSpike_Generator : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(TestSpike_Generator);
    CPPUNIT_TEST(testPeriodicSpikes);
//    CPPUNIT_TEST(testRandomSpikes);
//    CPPUNIT_TEST(testPoissonSpikes);
    CPPUNIT_TEST_SUITE_END();

    public:
        void setUp(void);
        void tearDown(void);

    protected:
        void testPeriodicSpikes(void);
//        void testRandomSpikes(void);
//        void testPoissonSpikes(void);

    private:
        Spike_Generator *mTestObj;
};

void TestSpike_Generator::testPeriodicSpikes(void)
{
    uint n = 10;
    uint step = 2;
    float value = 5;

    std::vector<float> vec = mTestObj->generate_periodic_spikes(n, step, value);
    std::vector<float> assert_vec = {0, value, 0, value, 0, value, 0, value, 0, value};

    for (uint i = 0; i < n; i++)
        CPPUNIT_ASSERT(vec[i] == assert_vec[i]);
    std::cout << "\n===== TestSpike_Generator : PERIODIC SPIKES VECTOR TEST PASSED =====" << std::endl;
}

void TestSpike_Generator::setUp(void)
{
    mTestObj = new Spike_Generator();
}

void TestSpike_Generator::tearDown(void)
{
    delete mTestObj;
}

CPPUNIT_TEST_SUITE_REGISTRATION( TestSpike_Generator );

int main(int argc, char* argv[])
{
    // informs test-listener about testresults
    CPPUNIT_NS::TestResult testresult;

    // register listener for collecting the test-results
    CPPUNIT_NS::TestResultCollector collectedresults;
    testresult.addListener (&collectedresults);

    // register listener for per-test progress output
    CPPUNIT_NS::BriefTestProgressListener progress;
    testresult.addListener (&progress);

    // insert test-suite at test-runner by registry
    CPPUNIT_NS::TestRunner testrunner;
    testrunner.addTest (CPPUNIT_NS::TestFactoryRegistry::getRegistry().makeTest ());
    testrunner.run(testresult);

    // output results in compiler-format
    CPPUNIT_NS::CompilerOutputter compileroutputter(&collectedresults, std::cerr);
    compileroutputter.write ();

    // return 0 if tests were successful
    return collectedresults.wasSuccessful() ? 0 : 1;
}
