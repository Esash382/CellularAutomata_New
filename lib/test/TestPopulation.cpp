#include <iostream>
#include <cmath>
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

#include <dlfcn.h>
#include <string.h>

#include "../external/population.hpp"

#include <iostream>

using namespace CppUnit;
using namespace std;

typedef unsigned int uint;
typedef Population* (*populationCreatorFunction)(const char* filepath);

class TestPopulation: public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(TestPopulation);
    CPPUNIT_TEST(testNetworkSize);
    CPPUNIT_TEST(testComputeWeights);
    CPPUNIT_TEST_SUITE_END();

    public:
        void setUp(void);
        void tearDown(void);

    protected:
        void testNetworkSize(void);
        void testComputeWeights(void);
        void testSynapticBlock(void);

    private:
        void set_weight_matrix(void);

    private:
        void* handle;
        Population *mTestObj;
};

void TestPopulation::testNetworkSize(void)
{
    auto matrix = mTestObj->get_weight_matrix();
    CPPUNIT_ASSERT(mTestObj->get_network_size() == matrix.size());
    CPPUNIT_ASSERT(mTestObj->get_network_size() == matrix[0].size());
    std::cout << "\n===== TestPopulation : CHECK NETWORK SIZE TEST PASSED =====" << std::endl;
}

void TestPopulation::set_weight_matrix(void)
{
    std::vector<std::vector<double>> matrix = {
                {0, 2, 15, 15 },
                {2, 0, 15, 15 },
                {2, 2, 0,  0  },
                {2, 2, 0,  0  },
           };
    mTestObj->set_weight_matrix(matrix);
}

void TestPopulation::testComputeWeights(void)
{
    // for random external input case
    // for periodic external input case
    // for poisson external input case
    // for random stimulation

    // In case of external inputs, no actual synapses are 'on'
    // The pseudo-synapses are 'on'.
    mTestObj->population_network[0]->external_input = RANDOM;
    mTestObj->population_network[1]->external_input = RANDOM;
    double sum = 0.0;
    if (mTestObj->population_network[0]->m_ntwk_name == "ex") {
        uint size = mTestObj->population_network[0]->m_N;
        uint i = 1;
        uint n = rand() % mTestObj->time_vec.size();

        mTestObj->m_n_matrix[0] = ON;
        mTestObj->m_nf_matrix[0] = mTestObj->time_vec[n];
        sum = mTestObj->call_compute_weights(mTestObj->population_network[0], n, i);

        // weight of projection: 2.0 + external input value: 0.0
        CPPUNIT_ASSERT(sum == 2.0);

        mTestObj->population_network[0]->external_input_vec[n] = 3.0;
        sum = mTestObj->call_compute_weights(mTestObj->population_network[0], n, i);
        CPPUNIT_ASSERT(sum == 5.0);
        CPPUNIT_ASSERT(sum > mTestObj->population_network[0]->threshold);
        mTestObj->m_n_matrix[1] = ON;
    }

    sum = 0.0;
    if (mTestObj->population_network[0]->m_ntwk_name == "in") {
        uint size = mTestObj->population_network[1]->m_N;
        uint i = 3;
        uint n = rand() % mTestObj->time_vec.size();

        mTestObj->m_n_matrix[2] = ON;
        mTestObj->m_nf_matrix[2] = mTestObj->time_vec[n];
        sum = mTestObj->call_compute_weights(mTestObj->population_network[1], n, i);

        // weight of projection: 2.0 + external input value: 0.0
        CPPUNIT_ASSERT(sum == 4.0);

        mTestObj->population_network[1]->external_input_vec[n] = 3.0;
        sum = mTestObj->call_compute_weights(mTestObj->population_network[1], n, i);
        CPPUNIT_ASSERT(sum == 7.0);
    }

    // RESET NETWORK
    uint n = rand() % mTestObj->time_vec.size();
    mTestObj->call_threshold_block(mTestObj->population_network[0], n, 0, OFF);
    mTestObj->call_threshold_block(mTestObj->population_network[0], n, 1, OFF);
    mTestObj->call_threshold_block(mTestObj->population_network[1], n, 0, OFF);
    mTestObj->m_nf_matrix[0] = 0.0;
    mTestObj->m_nf_matrix[1] = 0.0;
    mTestObj->m_nf_matrix[2] = 0.0;
    CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == OFF);
    CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == OFF);
    CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == OFF);

    std::cout << "\n===== TestPopulation : COMPUTE WEIGHTS TESTS PASSED =====" << std::endl;

    // In case of activating a single neuron:
    // Switch 'on' all synapses of that neuron
    n = rand() % mTestObj->time_vec.size();
    mTestObj->population_network[0]->tau_del = 0.0f;
    mTestObj->population_network[1]->tau_del = 0.0f;
    mTestObj->call_activate_single_neuron(mTestObj->population_network[0], n, 1);
    CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][1] == mTestObj->time_vec[n]);
    CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][1] == 0.0);
    CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][1] == mTestObj->time_vec[n]);
    CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][1] == mTestObj->time_vec[n]);

    std::cout << "\n===== TestPopulation : SINGLE NEURON ACTIVATION TEST PASSED =====" << std::endl;

    mTestObj->call_activate_synapses(n);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][0] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][0] == S_ON);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][0] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][0] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][1] == S_ON);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][1] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][1] == S_ON);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][1] == S_ON);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][2] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][2] == S_ON);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][2] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][2] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][3] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][3] == S_ON);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][3] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][3] == S_OFF);

    std::cout << "\n===== TestPopulation : ACTIVATE SYNAPSES WITH NO DELAY TEST PASSED =====" << std::endl;

    mTestObj->call_deactivate_synapses(n+100);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][0] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][0] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][0] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][0] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][1] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][1] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][1] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][1] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][2] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][2] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][2] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][2] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][3] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][3] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][3] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][3] == S_OFF);

    CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][0] == 0.0f);
    CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][0] == 0.0f);
    CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][0] == 0.0f);
    CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][0] == 0.0f);
    CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][1] == 0.0f);
    CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][1] == 0.0f);
    CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][1] == 0.0f);
    CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][1] == 0.0f);
    CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][2] == 0.0f);
    CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][2] == 0.0f);
    CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][2] == 0.0f);
    CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][2] == 0.0f);
    CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][3] == 0.0f);
    CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][3] == 0.0f);
    CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][3] == 0.0f);
    CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][3] == 0.0f);

    std::cout << "\n===== TestPopulation : DEACTIVATE SYNAPSES TEST PASSED =====" << std::endl;

    // Reset network properties
    mTestObj->population_network[0]->tau_del = 5.0f;
    mTestObj->population_network[1]->tau_del = 5.0f;

    n = rand() % mTestObj->time_vec.size();

    mTestObj->call_activate_single_neuron(mTestObj->population_network[0], n, 1);
    mTestObj->call_activate_synapses(n);

    CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][0] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][0] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][0] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][0] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][1] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][1] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][1] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][1] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][2] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][2] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][2] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][2] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][3] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][3] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][3] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][3] == S_OFF);

    mTestObj->call_activate_synapses(n+50);
    CPPUNIT_ASSERT(std::abs(mTestObj->m_sf_matrix[0][1] - mTestObj->time_vec[n+50]) < epsilon);
    CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][1] == 0.0);
    CPPUNIT_ASSERT(std::abs(mTestObj->m_sf_matrix[2][1] - mTestObj->time_vec[n+50]) < epsilon);
    CPPUNIT_ASSERT(std::abs(mTestObj->m_sf_matrix[3][1] == mTestObj->time_vec[n+50]) < epsilon);

    CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][0] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][0] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][0] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][0] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][1] == S_ON);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][1] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][1] == S_ON);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][1] == S_ON);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][2] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][2] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][2] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][2] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][3] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][3] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][3] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][3] == S_OFF);

    mTestObj->call_activate_synapses(n+100);

    CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][0] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][0] == S_ON);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][0] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][0] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][1] == S_ON);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][1] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][1] == S_ON);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][1] == S_ON);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][2] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][2] == S_ON);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][2] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][2] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][3] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][3] == S_ON);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][3] == S_OFF);
    CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][3] == S_OFF);

    std::cout << "\n===== TestPopulation : ACTIVATE SYNAPSES AFTER DELAY TEST PASSED =====" << std::endl;
}

void TestPopulation::setUp(void)
{
	handle = dlopen("libcaf.so", RTLD_NOW);
	char* error;
    const char* filepath = "/home/esash/Documents/Projects/CellularAutomata/CellularAutomata_Fast/";

	if (!handle) {
		fputs (dlerror(), stderr);
		exit(1);
    }

	populationCreatorFunction create = (populationCreatorFunction) dlsym(handle, "ext_create_population");
	if ((error = dlerror()) != NULL)  {
		fputs(error, stderr);
		exit(1);
	}

    mTestObj = (*create)(filepath);

	if (!mTestObj) {
		dlclose(handle);
		return;
	}

    mTestObj->add_network(filepath, "example");
    set_weight_matrix();
}

void TestPopulation::tearDown(void)
{
	free(mTestObj);
    dlclose(handle);
}

CPPUNIT_TEST_SUITE_REGISTRATION( TestPopulation );

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
