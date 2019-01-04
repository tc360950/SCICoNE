//
// Created by Tuncel  Mustafa Anil on 1/3/19.
//

#include <cxxopts.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "SignalProcessing.h"
#include "Utils.h"

using namespace std;

vector<double> breakpoint_detection(vector<vector<double>> &mat, int window_size = 5);

int main( int argc, char* argv[]) {

    int n_cells;
    int n_bins = 10000;
    int window_size = 10;
    double threshold_coefficient = 3.0;
    string f_name_postfix;
    string region_sizes_file;
    string d_matrix_file;

    cxxopts::Options options("Breakpoint detection executable", "detects the breakpoints in the genome across all cells.");
    options.add_options()
            ("d_matrix_file", "Path to the counts matrix file, delimiter: ' ', line separator: '\n' ", cxxopts::value(d_matrix_file))
            ("n_bins", "Number of bins in the input matrix", cxxopts::value(n_bins))
            ("n_cells", "Number of cells in the input matrix", cxxopts::value(n_cells))
            ("window_size", "the size of the window used in breakpoint detection", cxxopts::value(window_size))
            ("threshold", "the coefficient of the breakpoint threshold", cxxopts::value(threshold_coefficient))
            ("postfix", "Postfix to be added to the output files, this is useful when you are running multiple simulations through a work flow management system", cxxopts::value(f_name_postfix))
            ;

    auto result = options.parse(argc, argv);

    if (not result.count("d_matrix_file"))
    {
        cerr << "the D matrix file is not provided."<<endl;
        return EXIT_FAILURE;
    }
    if (not result.count("n_bins"))
    {
        cerr << "the number of bins is not provided."<<endl;
        return EXIT_FAILURE;
    }
    if (not result.count("n_cells"))
    {
        cerr << "the number of cells is not provided."<<endl;
        return EXIT_FAILURE;
    }

    // read the input d_matrix
    vector<vector<double>> d_bins(n_cells, vector<double>(n_bins));
    Utils::read_counts(d_bins, d_matrix_file);

    // create the region_sizes
    vector<int> region_sizes;
    vector<vector<double>> d_regions;

    // perform segmentation and peak detection, then define the region sizes
    SignalProcessing dsp;

    vector<double> s_p = breakpoint_detection(d_bins, window_size);

    vector<double> sp_cropped = dsp.crop(s_p, window_size);

    // median normalise sp_cropped
    dsp.median_normalise(sp_cropped);
    vector<double> sp_cropped_copy(sp_cropped); // the copy sp vector that'll contain the NaN values

    int lb = 0;
    size_t ub = sp_cropped.size();

    vector<int> all_max_ids;
    map<double, pair<unsigned,unsigned>> q_map; // a map that serves as a queue
    // use a map because you want the values to be always sorted
    // start with the 0.0 value, it will be removed from the map before other values are inserted
    q_map.emplace(0.0 , std::make_pair(lb,ub)); // initial boundaries

    while(!q_map.empty())
    {
        pair<unsigned,unsigned> elem = q_map.rbegin()->second; // use the rbegin to get the largest val
        q_map.erase(q_map.rbegin()->first); // remove the largest val by it's key

        if (elem.second - elem.first > window_size)
        {
            int max_idx = dsp.evaluate_peak(sp_cropped, sp_cropped_copy, elem.first, elem.second,
                                            threshold_coefficient);
            if (max_idx != -1) // -1 means rejected
            {

                // replace the nearby bins by nan
                int start_idx, stop_idx;
                start_idx = max_idx - window_size;
                stop_idx = max_idx + window_size;

                // check the boundries
                if (start_idx < 0)
                    start_idx = 0;
                if (stop_idx > sp_cropped.size())
                    stop_idx = sp_cropped.size();
                // set the nearby bins to nan
                for (int i = start_idx; i < stop_idx; ++i) {
                    sp_cropped_copy[i] = std::nan("");
                }

                // compute the left and right medians
                vector<double> left_vec(sp_cropped.begin() + elem.first, sp_cropped.begin() + max_idx);
                vector<double> right_vec(sp_cropped.begin() + max_idx + 1, sp_cropped.begin() + elem.second);

                if (!left_vec.empty())
                {
                    double median_left = MathOp::median(left_vec);
                    // normalise bins on the left by left median
                    for (int i = elem.first; i < max_idx; ++i) {
                        sp_cropped[i] /= median_left;
                    }
                    if (max_idx - elem.first > window_size)
                    {
                        try
                        {
                            int max_left_idx = dsp.find_highest_peak(sp_cropped, elem.first, max_idx);
                            q_map.emplace(sp_cropped[max_left_idx] , std::make_pair(elem.first,max_idx));
                        }
                        catch (const std::runtime_error& e)
                        {
                            std::cerr << e.what() << std::endl;
                        }
                    }
                }
                if (!right_vec.empty())
                {
                    double median_right = MathOp::median(right_vec);
                    // normalise bins on the right by right median
                    for (int j = max_idx + 1; j < elem.second; ++j) {
                        sp_cropped[j] /= median_right;
                    }
                    if (elem.second - (max_idx+1) > window_size)
                    {
                        try
                        {
                            int max_right_idx = dsp.find_highest_peak(sp_cropped, max_idx+1, elem.second);
                            q_map.emplace(sp_cropped[max_right_idx] , std::make_pair(max_idx + 1,elem.second));
                        }
                        catch (const std::runtime_error& e)
                        {
                            std::cerr << e.what() << std::endl;
                        }
                    }
                }

                all_max_ids.push_back(max_idx);

            }
        }
    }

    std::sort(all_max_ids.begin(), all_max_ids.end());

    std::ofstream tree_file("./"+ f_name_postfix+"_segmented_regions.txt");

    for (int k = 0; k < all_max_ids.size(); ++k) {
        // write it to file
        tree_file << all_max_ids[k] + window_size << endl;
    }

    std::ofstream reg_sizes_file("./"+ f_name_postfix+"_segmented_region_sizes.txt");
    int cum_sum = 0;

    for (int k = 0; k < all_max_ids.size(); ++k) {
        // write it to file
        if (k == 0)
            reg_sizes_file << all_max_ids[k] + window_size << endl;
        else
            reg_sizes_file << all_max_ids[k] - all_max_ids[k-1] << endl;
    }

    return EXIT_SUCCESS;
}

vector<double> breakpoint_detection(vector<vector<double>> &mat, int window_size)
{
    /*
     * Performs the breakpoint detection
     * */

    size_t n_cells = mat.size();

    // compute the AIC scores

    vector<vector<double>> aic_vec = MathOp::likelihood_ratio(mat,window_size);

    vector<vector<double>> sigma;

    size_t n_breakpoints = aic_vec.size();
    cout <<"n_breakpoints: " << n_breakpoints << " n_cells: " << n_cells <<endl;

    for (auto &vec: aic_vec) // compute sigma matrix
    {
        auto res = MathOp::combine_scores(vec);
        sigma.push_back(res);
    }

    vector<double> log_priors;
    for (int j = 0; j < n_cells; ++j) {
        log_priors.push_back(MathOp::breakpoint_log_prior(j, n_cells,0.5));
    }


    vector<vector<long double>> log_posterior;

    for (int k = 0; k < n_breakpoints; ++k) {
        log_posterior.push_back(vector<long double>());
        for (int j = 0; j < n_cells; ++j) {
            long double val = log_priors[j] + sigma[k][j];
            log_posterior[k].push_back(val);
        }
    }

    vector<vector<long double>> posterior;
    int k_star = 4;

    vector<double> s_p;

    for (int l = 0; l < n_breakpoints; ++l)
    {
        posterior.push_back(vector<long double>());

        long double max_num = *max_element(log_posterior[l].begin(), log_posterior[l].begin()+k_star-1);
        long double max_denom = *max_element(log_posterior[l].begin(), log_posterior[l].end());

        for (int j = 0; j < k_star - 1; ++j) {
            long double val =exp(log_posterior[l][j] - max_num);
            posterior[l].push_back(val);
        }
        for (int k = k_star -1 ; k < log_posterior[l].size(); ++k) {
            long double val =exp(log_posterior[l][k] - max_denom);
            posterior[l].push_back(val);
        }


        long double sp_num = std::accumulate(posterior[l].begin(), posterior[l].begin()+k_star-1, 0.0);
        sp_num  = log(sp_num) + max_num;
        long double sp_denom = std::accumulate(posterior[l].begin(), posterior[l].end(), 0.0);
        sp_denom = log(sp_denom) + max_denom;

        double sp_val = sp_denom-sp_num;

        s_p.push_back(sp_val);

    }

    return s_p;

}