"""
Process Rat Tumor trial data and fit Bayesian model to it. All to be done from command line.
Any plotting of results will be handled externally, probably in a Python-coded Jupyter
"""

using CSV
using Printf
using ArgParse
using DataFrames
using Statistics
using SpecialFunctions
using JSON

#############################

"""
Extract key arguments for the scipt
using https://argparsejl.readthedocs.io/en/latest/argparse.html

Returns:
--------
    Dict str->Any with command line args
"""
function parse_cmd()
    arg_parse_settings = ArgParseSettings()
    
    @add_arg_table! arg_parse_settings begin
        "--src_data_csv"
            help = "Path to CSV with rat tumor trials trials with columns tumor_count, rat_count"
            arg_type = String
            
         "--save_json"
            help = "Path to JSON file where all the outputs will be saved (data or pointers to data)"
            arg_type = String
            
         "--eta_min"
            help = "Minimum value of eta"
            arg_type = Float64
            
         "--eta_max"
            help = "Maximum value of eta"
            arg_type = Float64
            
         "--eta_step_count"
            help = "Number of steps along eta"
            arg_type = Int
            
         "--nu_min"
            help = "Minimum value of nu"
            arg_type = Float64
            
         "--nu_max"
            help = "Maximum value of nu"
            arg_type = Float64
            
         "--nu_step_count"
            help = "Number of steps along nu"
            arg_type = Int
    end
    
    return parse_args(arg_parse_settings)
end


#############################

"""
Given command, line arguments return array of scatter
values for eta and nu that cover the region, a grid that is

Aruguments:
-----------
    Dictionary with keys `nu_min`, `nu_max`, `nu_step_count` and and same for eta
    will be used to create 2 1d scatter arrays, which are flattened grid values for
    eta-nu grid
    
Returns:
--------
    nu_scatter_array, eta_scatter_array 1d arrays of the same size that cover the 
    grid specified by the input params
"""
function get_nu_eta_scatter_arr(
    cmd_args
)::Tuple{Array{Float64, 1}, Array{Float64, 1}}
    scatter_count = cmd_args["nu_step_count"] * cmd_args["eta_step_count"]
    
    # allocate memory
    nu_scat_arr = Array{Float64}(undef, scatter_count)
    eta_scat_arr = Array{Float64}(undef, scatter_count)
    
    # prepare indivudual range
    nu_rng = LinRange(cmd_args["nu_min"], cmd_args["nu_max"], cmd_args["nu_step_count"])
    eta_rng = LinRange(cmd_args["eta_min"], cmd_args["eta_max"], cmd_args["eta_step_count"])
    
    # loop through both ranges
    for (i_scat, (nu_val, eta_val))=enumerate(Iterators.product(nu_rng, eta_rng))
        nu_scat_arr[i_scat] = nu_val
        eta_scat_arr[i_scat] = eta_val
    end
    
    return nu_scat_arr, eta_scat_arr
end

#############################

"""
Following gelman compute joint posterior for alpha,beta whilst working on
nu=log(alpha / beta), eta=log(alpha+beta) grid

Arguments:
-----------
    nu_scat_arr: scatter point coordinates where eta-nu space has to be sampled
    eta_scat_arr: scatter point coordinates where eta-nu space has to be sampled
    rat_counts_arr: array of rat counts in each successive experiment
    tumor_counts arr: array of tumors in rats in each successive rat
    
Returns:
--------
    Logarithm of normalized (to sum=1) posterior probability, on eta-nu grid, sampled at points
    specified by nu_scar_arr and eta_scat_arr
"""
function log_posterior_on_nu_eta_grid(
    nu_scat_arr,
    eta_scat_arr,
    tumor_counts_arr,
    rat_counts_arr
)
    # compute alpha & beta for each point on the scatter
    alpha_scat_arr = exp.( nu_scat_arr .+ eta_scat_arr ) ./ ( 1. .+ exp.(nu_scat_arr) )
    beta_scat_arr = exp.(eta_scat_arr) ./ ( 1. .+ exp.(nu_scat_arr) )
    
    # add prior and normalize (to avoid overflow)
    log_post_scat_arr = (-5. / 2.) * log.(alpha_scat_arr .+ beta_scat_arr)
    log_post_scat_arr = log_post_scat_arr .+ log.(alpha_scat_arr )
    log_post_scat_arr = log_post_scat_arr .+ log.(beta_scat_arr )
    #
    norm_log_post(lp_arr) = lp_arr .- median(lp_arr)  
    log_post_scat_arr = norm_log_post(log_post_scat_arr)
    
    # now loop through the data points
    # adjusting posterior each time
    for (cur_tumor_count, cur_rat_count)=zip(tumor_counts_arr, rat_counts_arr)
        log_post_scat_arr = log_post_scat_arr .+ loggamma.(alpha_scat_arr .+ beta_scat_arr)
        log_post_scat_arr = log_post_scat_arr .- loggamma.(alpha_scat_arr)
        log_post_scat_arr = log_post_scat_arr .- loggamma.(beta_scat_arr)
        
        log_post_scat_arr = log_post_scat_arr .+ loggamma.(alpha_scat_arr .+ cur_tumor_count)
        log_post_scat_arr = log_post_scat_arr .+ loggamma.(beta_scat_arr .+ cur_rat_count .- cur_tumor_count)
        log_post_scat_arr = log_post_scat_arr .- loggamma.(alpha_scat_arr .+ beta_scat_arr .+ cur_rat_count)
        
        log_post_scat_arr = norm_log_post(log_post_scat_arr)
    end
    
    # normalize to 1
    cur_sum = sum(exp.(log_post_scat_arr))
    log_post_scat_arr = log_post_scat_arr .- log(cur_sum)
    
    return log_post_scat_arr
end

#############################

"""
Entry point to the program
"""
function main()
    # extract the path to the trial CSV
    cmd_args = parse_cmd()
    
    # load data as a dataframe
    trial_df = CSV.read(cmd_args["src_data_csv"], DataFrame)
    
    # prepare eta, nu values
    nu_scat_arr, eta_scat_arr = get_nu_eta_scatter_arr(cmd_args)
    
    # get the log of the posterior
    log_post_scat_arr = log_posterior_on_nu_eta_grid(
        nu_scat_arr,
        eta_scat_arr,
        Array{Float64}(trial_df[:, "tumor_count"]),
        Array{Float64}(trial_df[:, "rat_count"])
    )
    
    # convert to posterior
    post_scat_arr = exp.(log_post_scat_arr)
    
    # saving result into JSON
    open(cmd_args["save_json"], "w") do fh
        JSON.print(fh, Dict(
            "nu_scat_arr" => nu_scat_arr,
            "eta_scat_arr" => eta_scat_arr,
            "post_scat_arr" => post_scat_arr,
            "status" => "done"
        ))
    end
    
    @printf("Done\n")
end

#########

main()
