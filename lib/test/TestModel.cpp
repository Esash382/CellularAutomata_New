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

class TestModel: public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(TestModel);
    CPPUNIT_TEST(testPeriodicExternalInputWithoutDelay);
    CPPUNIT_TEST(testPeriodicExternalInputWithDelay);
    CPPUNIT_TEST_SUITE_END();

    public:
        void setUp(void);
        void tearDown(void);

    protected:
        void testPeriodicExternalInputWithoutDelay(void);
        void testPeriodicExternalInputWithDelay(void);

    private:
        void* handle;
        Population *mTestObj;
};

void TestModel::testPeriodicExternalInputWithoutDelay(void)
{
    // EINPUT stimulation : PERIODIC
    std::vector<uint> ntwk_id = { 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0};
    std::vector<uint> neuron_id = { 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0};

    mTestObj->population_network[0]->external_input = PERIODIC;
    mTestObj->population_network[1]->external_input = PERIODIC;

    mTestObj->time_vec = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    mTestObj->total_time = mTestObj->time_vec.size();
    mTestObj->time_step = 1.0f;

    mTestObj->population_network[0]->external_input_value = 3;
    mTestObj->population_network[0]->ext_step = 1.0f;
    mTestObj->population_network[0]->external_input_vec = {3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 0};
    mTestObj->population_network[0]->tau_ref = 2.0f;
    mTestObj->population_network[0]->tau_ap = 1.0f;
    mTestObj->population_network[0]->tau_dur = 2.0f;
    mTestObj->population_network[0]->tau_del = 0.0f;
    mTestObj->population_network[0]->threshold = 1.0f;

    mTestObj->population_network[1]->external_input_value = 3;
    mTestObj->population_network[1]->ext_step = 1.0f;
    mTestObj->population_network[1]->external_input_vec = {3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 0};
    mTestObj->population_network[1]->tau_ref = 2.0f;
    mTestObj->population_network[1]->tau_ap = 1.0f;
    mTestObj->population_network[1]->tau_dur = 2.0f;
    mTestObj->population_network[1]->tau_del = 0.0f;
    mTestObj->population_network[1]->threshold = 1.0f;

    std::vector<std::vector<double>> matrix = {
                {0, 2, 15, 15 },
                {2, 0, 15, 15 },
                {2, 2, 0,  0  },
                {2, 2, 0,  0  },
           };
    mTestObj->set_weight_matrix(matrix);

    // Check if all the states and firing times of synapses are 0
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

    srand(NULL);
    for (uint n = 0; n < mTestObj->time_vec.size(); n++) {
        mTestObj->call_deactivate_neurons(n);
        mTestObj->call_synaptic_block(mTestObj->population_network[ntwk_id[n]], n, neuron_id[n]);
        mTestObj->call_update_stats(n);

        // Check state matrix and synapse firing matrix for each neuron
        if ( n == 0 ) {
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][1] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][2] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][3] == S_ON);
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

            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == OFF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 0.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 1.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : network : " << mTestObj->population_network[ntwk_id[n]]->m_ntwk_name << " : neuron : " << neuron_id[n] << " : PASSED =====" << std::endl;  
        } else if ( n == 1 ) {
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][1] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][1] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][2] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][3] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][3] == S_OFF);

            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][0] == 1.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][1] == 1.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][3] == 0.0f);

            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == OFF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 1.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 0.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.5);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : network : " << mTestObj->population_network[ntwk_id[n]]->m_ntwk_name << " : neuron : " << neuron_id[n] << " : PASSED =====" << std::endl;  
        } else if ( n == 2 ) {
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][1] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][2] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][3] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][3] == S_OFF);

            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][0] == 2.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][0] == 1.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][1] == 1.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][2] == 2.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][3] == 2.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][3] == 0.0f);

            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == OFF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 1.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 0.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.5);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.5);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : network : " << mTestObj->population_network[ntwk_id[n]]->m_ntwk_name << " : neuron : " << neuron_id[n] << " : PASSED =====" << std::endl;  
        } else if ( n == 3 ) {
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][1] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][2] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][3] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][3] == S_OFF);

            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][0] == 2.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][0] == 3.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][1] == 3.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][2] == 2.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][3] == 2.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][3] == 0.0f);

            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == ON);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 1.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 3.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.5);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.5);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : network : " << mTestObj->population_network[ntwk_id[n]]->m_ntwk_name << " : neuron : " << neuron_id[n] << " : PASSED =====" << std::endl;  
        } else if ( n == 4 ) {
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][1] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][1] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][2] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][3] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][3] == S_OFF);

            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][0] == 3.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][1] == 4.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][1] == 3.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][2] == 4.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][3] == 4.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][3] == 0.0f);

            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 4.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 3.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.5);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.5);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : network : " << mTestObj->population_network[ntwk_id[n]]->m_ntwk_name << " : neuron : " << neuron_id[n] << " : PASSED =====" << std::endl;  
        } else if ( n == 5 ) {
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][1] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][1] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][2] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][3] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][3] == S_OFF);

            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][0] == 5.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][1] == 4.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][1] == 5.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][2] == 4.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][3] == 4.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][3] == 0.0f);

            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 4.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 5.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 3.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.5);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.5);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : network : " << mTestObj->population_network[ntwk_id[n]]->m_ntwk_name << " : neuron : " << neuron_id[n] << " : PASSED =====" << std::endl;  
        } else if ( n == 6 ) {
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][1] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][2] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][3] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][3] == S_OFF);

            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][0] == 6.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][0] == 5.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][1] == 5.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][2] == 6.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][3] == 6.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][3] == 0.0f);

            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == OFF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 4.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 6.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 5.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 0.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.5);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.5);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : network : " << mTestObj->population_network[ntwk_id[n]]->m_ntwk_name << " : neuron : " << neuron_id[n] << " : PASSED =====" << std::endl;  
        } else if ( n == 7 ) {
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][1] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][2] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][3] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][3] == S_OFF);

            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][0] == 6.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][0] == 7.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][1] == 7.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][2] == 6.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][3] == 6.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][3] == 0.0f);

            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == ON);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 6.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 5.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 7.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.5);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.5);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : network : " << mTestObj->population_network[ntwk_id[n]]->m_ntwk_name << " : neuron : " << neuron_id[n] << " : PASSED =====" << std::endl;  
        } else if ( n == 8 ) {
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][1] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][1] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][2] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][3] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][3] == S_OFF);

            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][0] == 7.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][1] == 8.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][1] == 7.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][2] == 8.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][3] == 8.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][3] == 0.0f);

            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 8.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 6.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 7.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.5);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.5);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : network : " << mTestObj->population_network[ntwk_id[n]]->m_ntwk_name << " : neuron : " << neuron_id[n] << " : PASSED =====" << std::endl;  
        } else if ( n == 9 ) {
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][1] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][1] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][2] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][3] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][3] == S_OFF);

            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][0] == 9.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][1] == 8.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][1] == 9.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][2] == 8.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][3] == 8.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][3] == 0.0f);

            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 8.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 9.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 7.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.5);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.5);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : network : " << mTestObj->population_network[ntwk_id[n]]->m_ntwk_name << " : neuron : " << neuron_id[n] << " : PASSED =====" << std::endl;  
        } else if ( n == 10 ) {
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][1] == S_ON);
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
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][0] == 9.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][1] == 9.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][3] == 0.0f);

            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == OFF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 8.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 9.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 0.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
           CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.5);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.5);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : network : " << mTestObj->population_network[ntwk_id[n]]->m_ntwk_name << " : neuron : " << neuron_id[n] << " : PASSED =====" << std::endl;  
        }
    }

    std::cout << "\n===== TestModel : PERIODIC EXTERNAL INPUT WITHOUT DELAY TEST PASSED =====" << std::endl;
}

void TestModel::testPeriodicExternalInputWithDelay(void)
{
    // EINPUT stimulation : RANDOM AND PERIODIC
    std::vector<uint> ntwk_id = { 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
    std::vector<uint> neuron_id = { 0, 0, 1, 1, 0, 0, 1, 1, 0, 0};

    mTestObj->population_network[0]->external_input = RANDOM;
    mTestObj->population_network[1]->external_input = PERIODIC;

    mTestObj->time_vec = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    mTestObj->total_time = mTestObj->time_vec.size();
    mTestObj->time_step = 1.0f;

    mTestObj->population_network[0]->external_input_value = 3;
    mTestObj->population_network[0]->ext_step = 1.0f;
    mTestObj->population_network[0]->external_input_vec = {5, 1, 2, 3, 4, 3, 2, 1, 0, 1};
    mTestObj->population_network[0]->tau_ref = 2.0f;
    mTestObj->population_network[0]->tau_ap = 1.0f;
    mTestObj->population_network[0]->tau_dur = 2.0f;
    mTestObj->population_network[0]->tau_del = 1.0f;
    mTestObj->population_network[0]->threshold = 1.0f;

    mTestObj->population_network[1]->external_input_value = 1;
    mTestObj->population_network[1]->ext_step = 1.0f;
    mTestObj->population_network[1]->external_input_vec = {3, 0, 3, 0, 3, 0, 3, 0, 3, 0};
    mTestObj->population_network[1]->tau_ref = 2.0f;
    mTestObj->population_network[1]->tau_ap = 1.0f;
    mTestObj->population_network[1]->tau_dur = 2.0f;
    mTestObj->population_network[1]->tau_del = 2.0f;
    mTestObj->population_network[1]->threshold = 1.0f;

    std::vector<std::vector<double>> matrix = {
                { 0,   2, 15, 15 },
                { 2,   0, 15, 15 },
                {-2, -2,  0,  0  },
                {-2, -2,  0,  0  },
           };
    mTestObj->set_weight_matrix(matrix);

    // Check if all the states and firing times of synapses are 0
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

    srand(NULL);
    for (uint n = 0; n < mTestObj->time_vec.size(); n++) {
        mTestObj->call_deactivate_neurons(n);
        mTestObj->call_synaptic_block(mTestObj->population_network[ntwk_id[n]], n, neuron_id[n]);
        mTestObj->call_update_stats(n);

        // Check state matrix and synapse firing matrix for each neuron
        if ( n == 0 ) {
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
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][1] == 1.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][2] == 1.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][3] == 1.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][3] == 0.0f);

            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == OFF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 0.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 1.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : network : " << mTestObj->population_network[ntwk_id[n]]->m_ntwk_name << " : neuron : " << neuron_id[n] << " : PASSED =====" << std::endl;  
        } else if ( n == 1 ) {
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][1] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][2] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][3] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][3] == S_OFF);

            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][0] == 2.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][0] == 3.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][0] == 3.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][1] == 1.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][1] == 3.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][1] == 3.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][2] == 1.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][2] == 2.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][3] == 1.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][3] == 2.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][3] == 0.0f);

            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == ON);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 1.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 1.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 1.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.5);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 1.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : network : " << mTestObj->population_network[ntwk_id[n]]->m_ntwk_name << " : neuron : " << neuron_id[n] << " : PASSED =====" << std::endl;  
        } else if ( n == 2 ) {
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][1] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][2] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][2] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][3] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][3] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][3] == S_OFF);

            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][0] == 2.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][0] == 3.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][0] == 3.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][1] == 1.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][1] == 3.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][1] == 3.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][2] == 1.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][2] == 2.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][3] == 1.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][3] == 2.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][3] == 0.0f);

            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 1.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 1.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 1.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 1.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 1.0);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : network : " << mTestObj->population_network[ntwk_id[n]]->m_ntwk_name << " : neuron : " << neuron_id[n] << " : PASSED =====" << std::endl;  
        } else if ( n == 3 ) {
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][1] == S_OFF);
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

            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][0] == 2.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][0] == 3.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][0] == 3.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][1] == 3.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][1] == 3.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][2] == 2.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][3] == 2.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][3] == 0.0f);

            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 1.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 1.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 1.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.5);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 1.0);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : network : " << mTestObj->population_network[ntwk_id[n]]->m_ntwk_name << " : neuron : " << neuron_id[n] << " : PASSED =====" << std::endl;  
        } else if ( n == 4 ) {
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][1] == S_OFF);
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

            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][0] == 3.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][0] == 3.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][1] == 3.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][1] == 3.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][3] == 0.0f);

            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == OFF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 0.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 1.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 1.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : network : " << mTestObj->population_network[ntwk_id[n]]->m_ntwk_name << " : neuron : " << neuron_id[n] << " : PASSED =====" << std::endl;  
        } else if ( n == 5 ) {
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

            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == OFF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 0.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 1.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 1.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : network : " << mTestObj->population_network[ntwk_id[n]]->m_ntwk_name << " : neuron : " << neuron_id[n] << " : PASSED =====" << std::endl;  
        } else if ( n == 6 ) {
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
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][0] == 7.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][2] == 7.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][3] == 7.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][3] == 0.0f);

            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == OFF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 6.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 0.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 1.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : network : " << mTestObj->population_network[ntwk_id[n]]->m_ntwk_name << " : neuron : " << neuron_id[n] << " : PASSED =====" << std::endl;  
        } else if ( n == 7 ) {
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][2] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][3] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][3] == S_OFF);

            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][0] == 7.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][0] == 9.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][0] == 9.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][1] == 8.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][1] == 9.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][1] == 9.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][2] == 8.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][2] == 7.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][3] == 8.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][3] == 7.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][3] == 0.0f);

            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == ON);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 7.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 6.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 7.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 7.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.5);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 1.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : network : " << mTestObj->population_network[ntwk_id[n]]->m_ntwk_name << " : neuron : " << neuron_id[n] << " : PASSED =====" << std::endl;  
        } else if ( n == 8 ) {
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][1] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][2] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][2] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][3] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][3] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][3] == S_OFF);

            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][0] == 7.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][0] == 9.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][0] == 9.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][1] == 8.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][1] == 9.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][1] == 9.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][2] == 8.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][2] == 7.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][3] == 8.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][3] == 7.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][3] == 0.0f);

            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 7.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 6.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 7.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 7.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 1.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 1.0);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : network : " << mTestObj->population_network[ntwk_id[n]]->m_ntwk_name << " : neuron : " << neuron_id[n] << " : PASSED =====" << std::endl;  
        } else if ( n == 9 ) {
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][0] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][0] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][1] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][1] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][1] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][1] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][2] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][2] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[0][3] == S_ON);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[1][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[2][3] == S_OFF);
            CPPUNIT_ASSERT(mTestObj->m_s_matrix[3][3] == S_OFF);

            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][0] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][0] == 9.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][0] == 9.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][1] == 8.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][1] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][1] == 9.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][1] == 9.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][2] == 8.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][2] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[0][3] == 8.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[1][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[2][3] == 0.0f);
            CPPUNIT_ASSERT(mTestObj->m_sf_matrix[3][3] == 0.0f);

            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 7.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 7.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 7.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.5);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 1.0);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : network : " << mTestObj->population_network[ntwk_id[n]]->m_ntwk_name << " : neuron : " << neuron_id[n] << " : PASSED =====" << std::endl;  
        }
    }

    std::cout << "\n===== TestModel : PERIODIC EXTERNAL INPUT WITH DELAY TEST PASSED =====" << std::endl;
}

void TestModel::setUp(void)
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
}

void TestModel::tearDown(void)
{
	free(mTestObj);
    dlclose(handle);
}

CPPUNIT_TEST_SUITE_REGISTRATION( TestModel );

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
