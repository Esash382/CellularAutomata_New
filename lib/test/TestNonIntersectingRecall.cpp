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

class TestNonIntersectionRecall: public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(TestNonIntersectionRecall);
    CPPUNIT_TEST(testRecallNetwork);
    CPPUNIT_TEST_SUITE_END();

    public:
        void setUp(void);
        void tearDown(void);

    protected:
        void testRecallNetwork(void);

    private:
        void* handle;
        Population *mTestObj;
};

void TestNonIntersectionRecall::testRecallNetwork(void)
{
    srand(0);

    // Set random number of external pseudo neurons to be activated at each time step
    {
        for (auto ntwk : mTestObj->population_network) {
            if (ntwk->external_input == RANDOM_TIME) {
                std::map<_time_t, std::vector<uint>> n_rand_activate;

                // which neurons
                std::vector<std::vector<unsigned int>> n_rand_indices = {
                                                            {9, 8, 2, 7, 5, 3, 1, 4, 0, 6},
                                                            {0, 3, 6, 7, 1, 9, 2, 5, 4, 8},
                                                            {0, 8, 9, 1, 4, 7, 6, 2, 5, 3},
                                                            {6, 5, 2, 3, 0, 1, 4, 8, 7, 9},
                                                            {1, 9, 2, 8, 4, 0, 5, 3, 7, 6},
                                                            {6, 4, 0, 9, 3, 5, 7, 8, 2, 1},
                                                            {6, 2, 3, 7, 1, 4, 8, 9, 0, 5},
                                                            {5, 7, 0, 1, 9, 8, 3, 4, 6, 2},
                                                            {7, 2, 6, 3, 1, 8, 4, 9, 5, 0},
                                                            {6, 5, 3, 2, 7, 4, 0, 8, 1, 9},
                                                            {6, 3, 9, 0, 7, 2, 4, 5, 8, 1},
                                                            {8, 4, 1, 0, 2, 9, 5, 7, 6, 3},
                                                            {0, 5, 2, 7, 1, 4, 8, 3, 6, 9},
                                                            {7, 0, 1, 9, 4, 5, 2, 8, 6, 3},
                                                            {3, 6, 0, 8, 4, 9, 5, 1, 2, 7},
                                                            {6, 4, 3, 8, 5, 2, 0, 7, 9, 1},
                                                            {1, 5, 2, 0, 3, 9, 4, 8, 7, 6},
                                                            {0, 5, 2, 3, 4, 9, 7, 8, 1, 6},
                                                            {9, 4, 3, 7, 5, 0, 8, 6, 1, 2},
                                                            {3, 9, 6, 5, 1, 7, 8, 0, 2, 4}
                };
                CPPUNIT_ASSERT(n_rand_indices.size() == mTestObj->time_vec.size());
                for (uint j = 0; j < mTestObj->time_vec.size(); j++) {
                    n_rand_activate[mTestObj->time_vec[j]] = n_rand_indices[j];
                }

                mTestObj->n_rand_map[ntwk->m_ntwk_id] = n_rand_activate;
            }
        }
    }

    std::vector<std::vector<double>> matrix = {
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 10, 10, 0, 10, 10, 0, 0}, 
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 10, 0, 10, 0, 10, 10, 0, 0, 0}, 
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 10, 0, 10, 10, 0, 10}, 
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 10, 0, 10, 0, 10, 0, 10}, 
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 10, 0, 10, 0, 10, 0, 10}, 
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 10, 10, 0, 10, 10}, 
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 10, 0, 0, 10, 0, 10, 10}, 
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 10, 0, 10, 0, 0, 10, 10, 0}, 
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 10, 0, 0, 10, 10, 10, 0, 0}, 
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 10, 10, 10, 0, 0, 10}, 
        {0, 50, 0, 50, 0, 50, 0, 50, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {50, 0, 50, 0, 50, 0, 50, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
        {0, 0, 50, 0, 50, 50, 0, 50, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
        {0, 0, 50, 0, 50, 0, 50, 0, 50, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
        {50, 50, 50, 0, 50, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
        {0, 50, 50, 0, 0, 50, 50, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
        {50, 0, 0, 0, 50, 0, 0, 50, 50, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
        {0, 0, 50, 0, 50, 50, 50, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
        {50, 50, 0, 50, 0, 0, 0, 50, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
        {50, 50, 50, 0, 0, 0, 50, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
        {-5, 0, -5, -5, 0, 0, -5, 0, 0, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
        {0, 0, 0, -5, 0, -5, -5, -5, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
        {0, -5, 0, -5, -5, 0, 0, 0, -5, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
        {-5, 0, -5, 0, -5, 0, 0, 0, -5, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
        {0, -5, -5, -5, -5, 0, 0, 0, 0, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
        {0, -5, -5, -5, -5, 0, 0, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
        {-5, -5, 0, -5, -5, 0, 0, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
        {0, -5, -5, -5, 0, 0, 0, -5, 0, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
        {0, -5, -5, -5, 0, -5, 0, 0, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
        {-5, 0, 0, -5, 0, -5, 0, -5, 0, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
    };
 
    mTestObj->set_weight_matrix(matrix);
    mTestObj->del_step = 0.0;

    for (uint n = 0; n < mTestObj->time_vec.size(); n++) {
        mTestObj->call_deactivate_neurons(n);

        bool flag = false;
        for (auto ntwk : mTestObj->population_network) {
            auto it = mTestObj->n_rand_map.find(ntwk->m_ntwk_id);
            if (it != mTestObj->n_rand_map.end()) {
                auto n_rand_activate = it->second;
                auto it2 = n_rand_activate.find(mTestObj->time_vec[n]);
                if (it2 != n_rand_activate.end()) {
                    auto neuron_vec = it2->second;
                    for (uint j = 0; j < neuron_vec.size(); j++) {
                       mTestObj->call_synaptic_block(ntwk, n, neuron_vec[j]);
                       flag = true;
                    }
                }
            }
        }

        if (!flag) {
            mTestObj->call_deactivate_synapses(n);
            mTestObj->call_activate_synapses(n);
            flag = false;
        }

        mTestObj->call_update_stats(mTestObj->time_vec[n]);
        // Check state matrix and synapse firing matrix for each neuron
        if ( n == 0 ) {
            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[4] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[5] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[6] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[7] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[8] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[9] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[10] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[11] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[12] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[13] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[14] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[15] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[16] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[17] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[18] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[19] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[20] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[21] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[22] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[23] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[24] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[25] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[26] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[27] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[28] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[29] == ON);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[4] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[5] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[6] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[7] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[8] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[9] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[10] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[11] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[12] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[13] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[14] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[15] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[16] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[17] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[18] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[19] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[20] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[21] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[22] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[23] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[24] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[25] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[26] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[27] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[28] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[29] == 0.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.9);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.1);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.0);

            // Check 'ext' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.8);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.2);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ac_stats[n] == 0.8);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->in_stats[n] == 0.2);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ref_stats[n] == 0.0);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : PASSED =====" << std::endl;  
        } else if ( n == 1 ) {
            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[4] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[5] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[6] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[7] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[8] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[9] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[10] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[11] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[12] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[13] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[14] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[15] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[16] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[17] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[18] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[19] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[20] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[21] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[22] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[23] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[24] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[25] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[26] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[27] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[28] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[29] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[4] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[5] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[6] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[7] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[8] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[9] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[10] == 1.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[11] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[12] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[13] == 1.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[14] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[15] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[16] == 1.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[17] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[18] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[19] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[20] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[21] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[22] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[23] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[24] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[25] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[26] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[27] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[28] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[29] == 0.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.1);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.9);

            // Check 'ext' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.3);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.7);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->in_stats[n] == 0.2);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ref_stats[n] == 0.8);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : PASSED =====" << std::endl;  
        } else if ( n == 2 ) {
            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[4] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[5] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[6] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[7] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[8] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[9] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[10] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[11] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[12] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[13] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[14] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[15] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[16] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[17] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[18] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[19] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[20] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[21] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[22] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[23] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[24] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[25] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[26] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[27] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[28] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[29] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[4] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[5] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[6] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[7] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[8] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[9] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[10] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[11] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[12] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[13] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[14] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[15] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[16] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[17] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[18] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[19] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[20] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[21] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[22] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[23] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[24] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[25] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[26] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[27] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[28] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[29] == 0.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.1);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.9);

            // Check 'ext' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.2);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.8);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ac_stats[n] == 0.1);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->in_stats[n] == 0.1);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ref_stats[n] == 0.8);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : PASSED =====" << std::endl;  
        } else if ( n == 3 ) {
            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[4] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[5] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[6] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[7] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[8] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[9] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[10] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[11] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[12] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[13] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[14] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[15] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[16] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[17] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[18] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[19] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[20] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[21] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[22] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[23] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[24] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[25] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[26] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[27] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[28] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[29] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[4] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[5] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[6] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[7] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[8] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[9] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[10] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[11] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[12] == 3.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[13] == 3.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[14] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[15] == 3.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[16] == 3.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[17] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[18] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[19] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[20] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[21] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[22] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[23] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[24] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[25] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[26] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[27] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[28] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[29] == 0.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 1.0);

            // Check 'ext' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.4);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.6);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->in_stats[n] == 0.1);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ref_stats[n] == 0.9);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : PASSED =====" << std::endl;  
        } else if ( n == 4 ) {
            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[4] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[5] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[6] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[7] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[8] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[9] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[10] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[11] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[12] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[13] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[14] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[15] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[16] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[17] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[18] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[19] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[20] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[21] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[22] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[23] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[24] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[25] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[26] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[27] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[28] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[29] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[4] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[5] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[6] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[7] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[8] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[9] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[10] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[11] == 4.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[12] == 4.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[13] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[14] == 4.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[15] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[16] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[17] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[18] == 4.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[19] == 4.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[20] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[21] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[22] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[23] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[24] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[25] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[26] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[27] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[28] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[29] == 0.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 1.0);

            // Check 'ext' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->in_stats[n] == 0.1);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ref_stats[n] == 0.9);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : PASSED =====" << std::endl;  
        } else if ( n == 5 ) {
            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[4] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[5] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[6] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[7] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[8] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[9] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[10] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[11] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[12] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[13] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[14] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[15] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[16] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[17] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[18] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[19] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[20] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[21] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[22] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[23] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[24] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[25] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[26] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[27] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[28] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[29] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[4] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[5] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[6] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[7] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[8] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[9] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[10] == 5.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[11] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[12] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[13] == 5.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[14] == 5.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[15] == 5.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[16] == 5.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[17] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[18] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[19] == 5.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[20] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[21] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[22] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[23] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[24] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[25] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[26] == 5.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[27] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[28] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[29] == 0.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 1.0);

            // Check 'ext' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.6);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.4);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ac_stats[n] == 0.1);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ref_stats[n] == 0.9);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : PASSED =====" << std::endl;  
        } else if ( n == 6 ) {
            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[4] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[5] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[6] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[7] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[8] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[9] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[10] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[11] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[12] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[13] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[14] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[15] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[16] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[17] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[18] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[19] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[20] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[21] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[22] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[23] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[24] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[25] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[26] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[27] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[28] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[29] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[4] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[5] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[6] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[7] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[8] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[9] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[10] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[11] == 6.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[12] == 6.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[13] == 6.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[14] == 6.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[15] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[16] == 6.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[17] == 6.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[18] == 6.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[19] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[20] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[21] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[22] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[23] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[24] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[25] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[26] == 5.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[27] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[28] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[29] == 0.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 1.0);

            // Check 'ext' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.7);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.3);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ref_stats[n] == 1.0);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : PASSED =====" << std::endl;  
        } else if ( n == 7 ) {
            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[4] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[5] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[6] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[7] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[8] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[9] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[10] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[11] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[12] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[13] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[14] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[15] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[16] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[17] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[18] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[19] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[20] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[21] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[22] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[23] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[24] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[25] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[26] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[27] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[28] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[29] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[4] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[5] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[6] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[7] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[8] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[9] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[10] == 7.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[11] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[12] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[13] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[14] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[15] == 7.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[16] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[17] == 7.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[18] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[19] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[20] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[21] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[22] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[23] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[24] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[25] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[26] == 5.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[27] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[28] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[29] == 0.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 1.0);

            // Check 'ext' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.3);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.7);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ref_stats[n] == 1.0);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : PASSED =====" << std::endl;  
        } else if ( n == 8 ) {
            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[4] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[5] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[6] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[7] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[8] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[9] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[10] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[11] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[12] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[13] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[14] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[15] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[16] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[17] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[18] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[19] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[20] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[21] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[22] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[23] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[24] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[25] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[26] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[27] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[28] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[29] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[4] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[5] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[6] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[7] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[8] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[9] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[10] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[11] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[12] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[13] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[14] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[15] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[16] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[17] == 8.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[18] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[19] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[20] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[21] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[22] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[23] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[24] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[25] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[26] == 5.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[27] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[28] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[29] == 0.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 1.0);

            // Check 'ext' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.1);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.9);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ref_stats[n] == 1.0);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : PASSED =====" << std::endl;  
        } else if ( n == 9 ) {
            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[4] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[5] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[6] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[7] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[8] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[9] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[10] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[11] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[12] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[13] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[14] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[15] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[16] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[17] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[18] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[19] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[20] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[21] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[22] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[23] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[24] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[25] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[26] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[27] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[28] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[29] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[4] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[5] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[6] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[7] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[8] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[9] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[10] == 9.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[11] == 9.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[12] == 9.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[13] == 9.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[14] == 9.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[15] == 9.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[16] == 9.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[17] == 9.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[18] == 9.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[19] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[20] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[21] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[22] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[23] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[24] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[25] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[26] == 5.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[27] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[28] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[29] == 0.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 1.0);

            // Check 'ext' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.9);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.1);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ref_stats[n] == 1.0);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : PASSED =====" << std::endl;  
        } else if ( n == 10 ) {
            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[4] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[5] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[6] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[7] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[8] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[9] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[10] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[11] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[12] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[13] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[14] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[15] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[16] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[17] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[18] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[19] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[20] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[21] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[22] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[23] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[24] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[25] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[26] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[27] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[28] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[29] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[4] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[5] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[6] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[7] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[8] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[9] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[10] == 10.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[11] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[12] == 10.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[13] == 10.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[14] == 10.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[15] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[16] == 10.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[17] == 10.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[18] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[19] == 10.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[20] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[21] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[22] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[23] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[24] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[25] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[26] == 5.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[27] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[28] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[29] == 0.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 1.0);

            // Check 'ext' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.7);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.3);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ref_stats[n] == 1.0);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : PASSED =====" << std::endl;  
        } else if ( n == 11 ) {
            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[4] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[5] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[6] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[7] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[8] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[9] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[10] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[11] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[12] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[13] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[14] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[15] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[16] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[17] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[18] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[19] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[20] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[21] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[22] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[23] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[24] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[25] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[26] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[27] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[28] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[29] == OFF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[4] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[5] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[6] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[7] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[8] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[9] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[10] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[11] == 11.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[12] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[13] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[14] == 11.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[15] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[16] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[17] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[18] == 11.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[19] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[20] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[21] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[22] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[23] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[24] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[25] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[26] == 5.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[27] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[28] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[29] == 0.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 1.0);

            // Check 'ext' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.3);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.7);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->in_stats[n] == 0.8);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ref_stats[n] == 0.2);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : PASSED =====" << std::endl;  
        } else if ( n == 12 ) {
            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[4] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[5] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[6] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[7] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[8] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[9] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[10] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[11] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[12] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[13] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[14] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[15] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[16] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[17] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[18] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[19] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[20] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[21] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[22] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[23] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[24] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[25] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[26] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[27] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[28] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[29] == OFF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[4] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[5] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[6] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[7] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[8] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[9] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[10] == 12.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[11] == 12.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[12] == 12.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[13] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[14] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[15] == 12.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[16] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[17] == 12.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[18] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[19] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[20] == 2.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[21] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[22] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[23] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[24] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[25] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[26] == 5.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[27] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[28] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[29] == 0.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.1);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.9);

            // Check 'ext' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->in_stats[n] == 0.8);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ref_stats[n] == 0.2);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : PASSED =====" << std::endl;  
        } else if ( n == 13 ) {
            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[4] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[5] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[6] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[7] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[8] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[9] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[10] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[11] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[12] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[13] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[14] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[15] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[16] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[17] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[18] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[19] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[20] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[21] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[22] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[23] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[24] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[25] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[26] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[27] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[28] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[29] == ON);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[4] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[5] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[6] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[7] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[8] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[9] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[10] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[11] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[12] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[13] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[14] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[15] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[16] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[17] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[18] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[19] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[20] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[21] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[22] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[23] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[24] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[25] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[26] == 5.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[27] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[28] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[29] == 13.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.1);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.9);

            // Check 'ext' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.8);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.2);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ac_stats[n] == 0.4);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->in_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ref_stats[n] == 0.1);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : PASSED =====" << std::endl;  
        } else if ( n == 14 ) {
            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[4] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[5] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[6] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[7] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[8] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[9] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[10] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[11] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[12] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[13] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[14] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[15] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[16] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[17] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[18] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[19] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[20] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[21] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[22] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[23] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[24] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[25] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[26] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[27] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[28] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[29] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[4] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[5] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[6] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[7] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[8] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[9] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[10] == 14.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[11] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[12] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[13] == 14.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[14] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[15] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[16] == 14.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[17] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[18] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[19] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[20] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[21] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[22] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[23] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[24] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[25] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[26] == 5.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[27] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[28] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[29] == 13.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 1.0);

            // Check 'ext' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.3);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.7);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->in_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ref_stats[n] == 0.5);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : PASSED =====" << std::endl;  
        } else if ( n == 15 ) {
            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[4] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[5] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[6] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[7] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[8] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[9] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[10] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[11] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[12] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[13] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[14] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[15] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[16] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[17] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[18] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[19] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[20] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[21] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[22] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[23] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[24] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[25] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[26] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[27] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[28] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[29] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[4] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[5] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[6] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[7] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[8] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[9] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[10] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[11] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[12] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[13] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[14] == 15.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[15] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[16] == 15.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[17] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[18] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[19] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[20] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[21] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[22] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[23] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[24] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[25] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[26] == 5.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[27] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[28] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[29] == 13.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 1.0);

            // Check 'ext' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.2);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.8);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->in_stats[n] == 0.5);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ref_stats[n] == 0.5);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : PASSED =====" << std::endl;  
        } else if ( n == 16 ) {
            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[4] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[5] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[6] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[7] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[8] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[9] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[10] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[11] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[12] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[13] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[14] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[15] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[16] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[17] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[18] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[19] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[20] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[21] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[22] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[23] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[24] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[25] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[26] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[27] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[28] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[29] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[4] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[5] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[6] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[7] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[8] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[9] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[10] == 16.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[11] == 16.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[12] == 16.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[13] == 16.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[14] == 16.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[15] == 16.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[16] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[17] == 16.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[18] == 16.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[19] == 16.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[20] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[21] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[22] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[23] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[24] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[25] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[26] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[27] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[28] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[29] == 13.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 1.0);

            // Check 'ext' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.9);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.1);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->in_stats[n] == 0.6);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ref_stats[n] == 0.4);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : PASSED =====" << std::endl;  
        } else if ( n == 17 ) {
            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[4] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[5] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[6] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[7] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[8] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[9] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[10] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[11] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[12] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[13] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[14] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[15] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[16] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[17] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[18] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[19] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[20] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[21] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[22] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[23] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[24] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[25] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[26] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[27] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[28] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[29] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[4] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[5] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[6] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[7] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[8] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[9] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[10] == 17.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[11] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[12] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[13] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[14] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[15] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[16] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[17] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[18] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[19] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[20] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[21] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[22] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[23] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[24] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[25] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[26] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[27] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[28] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[29] == 13.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.9);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.1);

            // Check 'ext' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.1);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.9);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ac_stats[n] == 0.0);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->in_stats[n] == 0.6);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ref_stats[n] == 0.4);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : PASSED =====" << std::endl;  
        } else if ( n == 18 ) {
            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[4] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[5] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[6] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[7] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[8] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[9] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[10] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[11] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[12] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[13] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[14] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[15] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[16] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[17] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[18] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[19] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[20] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[21] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[22] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[23] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[24] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[25] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[26] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[27] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[28] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[29] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 18.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 18.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 18.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[4] == 18.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[5] == 18.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[6] == 18.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[7] == 18.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[8] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[9] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[10] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[11] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[12] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[13] == 18.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[14] == 18.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[15] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[16] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[17] == 18.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[18] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[19] == 18.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[20] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[21] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[22] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[23] == 18.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[24] == 18.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[25] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[26] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[27] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[28] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[29] == 13.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.7);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.2);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.1);

            // Check 'ext' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.4);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.6);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ac_stats[n] == 0.2);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->in_stats[n] == 0.4);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ref_stats[n] == 0.4);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : PASSED =====" << std::endl;  
        } else if ( n == 19 ) {
            // States of neurons
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[0] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[1] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[2] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[3] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[4] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[5] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[6] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[7] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[8] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[9] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[10] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[11] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[12] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[13] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[14] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[15] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[16] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[17] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[18] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[19] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[20] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[21] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[22] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[23] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[24] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[25] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[26] == ON);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[27] == REF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[28] == OFF);
            CPPUNIT_ASSERT(mTestObj->m_n_matrix[29] == REF);

            // Firing times of neurons
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[0] == 18.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[1] == 18.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[2] == 18.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[3] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[4] == 18.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[5] == 18.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[6] == 18.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[7] == 18.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[8] == 19.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[9] == 19.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[10] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[11] == 19.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[12] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[13] == 19.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[14] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[15] == 19.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[16] == 19.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[17] == 19.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[18] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[19] == 19.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[20] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[21] == 19.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[22] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[23] == 18.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[24] == 18.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[25] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[26] == 19.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[27] == 13.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[28] == 0.0);
            CPPUNIT_ASSERT(mTestObj->m_nf_matrix[29] == 13.0);

            // Check 'ex' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ac_stats[n] == 0.2);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->in_stats[n] == 0.1);
            CPPUNIT_ASSERT(mTestObj->population_network[0]->ref_stats[n] == 0.7);

            // Check 'ext' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ac_stats[n] == 0.6);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->in_stats[n] == 0.4);
            CPPUNIT_ASSERT(mTestObj->population_network[1]->ref_stats[n] == 0.0);

            // Check 'in' network stats
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ac_stats[n] == 0.2);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->in_stats[n] == 0.2);
            CPPUNIT_ASSERT(mTestObj->population_network[2]->ref_stats[n] == 0.6);

            std::cout << "\n===== TestModel : " << __FUNCTION__ << " : time : " << n << " : PASSED =====" << std::endl;  
        }
    }

    std::cout << "\n===== TestNonIntersectionRecall : RANDOM NETWORKS TEST PASSED =====" << std::endl;
}

void TestNonIntersectionRecall::setUp(void)
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

    mTestObj->add_network(filepath, "example_recall_test");
}

void TestNonIntersectionRecall::tearDown(void)
{
	free(mTestObj);
    dlclose(handle);
}

CPPUNIT_TEST_SUITE_REGISTRATION( TestNonIntersectionRecall );

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
