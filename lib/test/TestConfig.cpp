#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <memory>
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

#include <string.h>

#include "../include/test.hpp"
#include "../include/config.hpp"

using namespace CppUnit;
using namespace std;

class TestConfig: public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(TestConfig);
    CPPUNIT_TEST(testGetFiles);
    CPPUNIT_TEST(testReadNeuronFile);
    CPPUNIT_TEST(testReadNetworkFile);
    CPPUNIT_TEST(testCreateNetwork);
    CPPUNIT_TEST_SUITE_END();

    public:
        void setUp(void);
        void tearDown(void);

    protected:
        void testGetFiles(void);
        void testReadNeuronFile(void);
        void testReadNetworkFile(void);
        void testCreateNetwork();

    private:
        Config *mTestObj;
};

void TestConfig::testGetFiles(void)
{
    std::vector<std::vector<string>> files = mTestObj->get_files();
    CPPUNIT_ASSERT(files.size() == 1);
    CPPUNIT_ASSERT(files[0].size() == 2);
    CPPUNIT_ASSERT(strncmp(files[0][0].c_str(), std::string("ex").c_str(), sizeof("ex")) == 0);
    CPPUNIT_ASSERT(strncmp(files[0][1].c_str(), std::string("in").c_str(), sizeof("in")) == 0);
    std::cout << "\n===== TestConfig : GET NETWORK FILES TEST PASSED =====" << std::endl;
}

void TestConfig::testReadNeuronFile(void)
{
    std::map<string, float> con = mTestObj->read_file("ex");
    CPPUNIT_ASSERT(con["id"] == 101.0);
    CPPUNIT_ASSERT(con["type"] == 0.0);
    CPPUNIT_ASSERT(con["N"] == 2.0);
    CPPUNIT_ASSERT(con["k"] == 2.0);
    CPPUNIT_ASSERT(con["threshold"] == 2.0);
    CPPUNIT_ASSERT(con["z"] == 1.0);
    CPPUNIT_ASSERT(con["ap"] == 1);
    CPPUNIT_ASSERT(con["ref"] == 10.0);
    CPPUNIT_ASSERT(con["del"] == 5.0);
    CPPUNIT_ASSERT(con["dur"] == 10.0);
    CPPUNIT_ASSERT(con["osc"] == 0.0);
    CPPUNIT_ASSERT(con["ext_type"] == 0);
    CPPUNIT_ASSERT(con["ext_val"] == 3);
    CPPUNIT_ASSERT(con["ext_step"] == 100);
    std::cout << "\n===== TestConfig : READ NEURON FILES TEST PASSED =====" << std::endl;
}

void TestConfig::testReadNetworkFile(void)
{
    std::map<string, float> con = mTestObj->read_file("ex_in");
    CPPUNIT_ASSERT(con["id"] == 201.0);
    CPPUNIT_ASSERT(con["ksd"] == 15.0);
    CPPUNIT_ASSERT(con["zsd"] == 2.0);
    CPPUNIT_ASSERT(con["kds"] == -1.0);
    CPPUNIT_ASSERT(con["zds"] == 2.0);
    std::cout << "\n===== TestConfig : READ NETWORK FILES TEST PASSED =====" << std::endl;
}

void TestConfig::testCreateNetwork(void)
{
    std::vector<std::vector<string>> files = mTestObj->get_files();
    shared_ptr<Network> nn = mTestObj->create_network(files[0][0]);
    CPPUNIT_ASSERT(nn->m_N == 2);
    std::cout << "\n===== TestConfig : CREATE NETWORK FILES TEST PASSED =====" << std::endl;
}

void TestConfig::setUp(void)
{
    mTestObj = new Config("/home/esash/Documents/Projects/CellularAutomata/CellularAutomata_Fast/", "example");
}

void TestConfig::tearDown(void)
{
    delete mTestObj;
}

CPPUNIT_TEST_SUITE_REGISTRATION( TestConfig );

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
