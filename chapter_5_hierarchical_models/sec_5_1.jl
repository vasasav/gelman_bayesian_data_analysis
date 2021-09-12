using CSV
using Printf
using ArgParse
using DataFrames
using Statistics

#############################

"""
Extract key arguments for the scipt
using https://argparsejl.readthedocs.io/en/latest/argparse.html
"""
function parse_cmd()
    arg_parse_settings = ArgParseSettings()
    
    @add_arg_table! arg_parse_settings begin
        "--trial_csv"
            help = "Path to CSV with trials with columns tumor_count, rat_count"
            arg_type = String
            default = "table_5_1_rat_tumors.csv"
    end
    
    return parse_args(arg_parse_settings)
end


"""
Print statistical summary of the data
"""
function print_trial_summary(trial_df)
    ### give some results on averages
    # 
    @printf(
        "\nAverage number of rats: %.3f\nAverage number of tumors: %.3f\n\n", 
        mean(trial_df.rat_count), 
        mean(trial_df.tumor_count)
    )
    
    
    ### key statistics
    mean_chance = mean(trial_df.tumor_count ./ trial_df.rat_count)
    var_chance = var(trial_df.tumor_count ./ trial_df.rat_count)
    
    @printf(
        "\nMean tumor chance: %.3e\nVariance: %.3e\n", 
        mean_chance,
        var_chance
    )
    
    ### work out the estimates for beta distribution
    beta_est = (mean_chance * ((1-mean_chance)^2) )/var_chance - (1 - mean_chance)
    alpha_est = mean_chance * beta_est / (1 - mean_chance)
    
    @printf("\nEquivalent beta-distribution params:\nalpha=%.3f\nbeta=%.3f\n\n", alpha_est, beta_est)
end


#############################

"""
Entry point to the program
"""
function main()
    # extract the path to the trial CSV
    cmd_args = parse_cmd()
    
    # read it as a dataframe
    trial_df = CSV.read(cmd_args["trial_csv"], DataFrame)

    # print some summary
    print_trial_summary(trial_df)
end

#########

main()
