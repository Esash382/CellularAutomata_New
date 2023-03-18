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
#include "../../include/plot.hpp"

#include <iostream>

using namespace CppUnit;
using namespace std;

typedef unsigned int uint;
typedef Population* (*populationCreatorFunction)(const char* filepath);

class TestRecallWeight: public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(TestRecallWeight);
    CPPUNIT_TEST(testRecallWeights);
    CPPUNIT_TEST_SUITE_END();

    public:
        void setUp(void);
        void tearDown(void);

    protected:
        void testRecallWeights(void);

    private:
        void set_p_rand_neuron_weight_matrix(void);

    private:
        void* handle;
        Population *mTestObj;
};

void TestRecallWeight::set_p_rand_neuron_weight_matrix(void)
{
    for (auto src_ntwk : mTestObj->population_network) {
        for (auto dst_ntwk : mTestObj->population_network) {
            if ((src_ntwk->m_type == PSEUDO_NEURON
                    && dst_ntwk->m_type != PSEUDO_NEURON)
                    && (src_ntwk->enable_learning == 1.0 
                        && dst_ntwk->enable_learning == 1.0)) {
                for (uint i = src_ntwk->start_from_row_index;
                          i < src_ntwk->start_from_row_index + src_ntwk->m_N;
                          i++) {
                    for (uint j = dst_ntwk->start_from_row_index;
                              j < dst_ntwk->start_from_row_index + dst_ntwk->m_N;
                              j++) {
                        if (mTestObj->m_w_matrix[i][j] != 0
                            && mTestObj->is_neuron_in_p_rand(i, src_ntwk->m_ntwk_id)
                            && mTestObj->is_neuron_in_p_rand(j, dst_ntwk->m_ntwk_id)) {
                            mTestObj->m_w_matrix[i][j] += dst_ntwk->learning_rate;
                        } else if (mTestObj->m_w_matrix[i][j] > 0) {
                            mTestObj->m_w_matrix[i][j] -= dst_ntwk->unlearning_rate;
//                            mTestObj->m_w_matrix[i][j] = 0;
                        }
                    }
                }
            }
        }
    }
}

void TestRecallWeight::testRecallWeights(void)
{
//    std::vector<std::vector<double>> matrix = {
//                {0, 2, 15, 15 },
//                {2, 0, 15, 15 },
//                {2, 2, 0,  0  },
//                {2, 2, 0,  0  },
//           };
//    mTestObj->set_weight_matrix(matrix);

    mTestObj->process_networks();
	Plot::plot_neurons_stats();
}

void TestRecallWeight::setUp(void)
{
	handle = dlopen("libca.so", RTLD_NOW);
	char* error;
    const char* filepath = "/home/ashraya/Documents/Notes/CellularAutomata_Fast";

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

    mTestObj->add_network(filepath, "example_recall");
    set_p_rand_neuron_weight_matrix();
}

void TestRecallWeight::tearDown(void)
{
	free(mTestObj);
    dlclose(handle);
}

CPPUNIT_TEST_SUITE_REGISTRATION( TestRecallWeight);

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
