/////////////////////////////////////////////////////////////////
// Contrafold.cpp
/////////////////////////////////////////////////////////////////

// include files
#ifdef MULTI
#include <mpi.h>
#endif
#include "Config.hpp"
#include "Options.hpp"
#include "Utilities.ipp"
#include "Utilities.hpp"
#include "ComputationWrapper.hpp"
#include "FileDescription.hpp"
#include "InferenceEngine.hpp"
#include "ParameterManager.hpp"
#include "OptimizationWrapper.hpp"

// constants
const double GAMMA_DEFAULT = 6;
const double KAPPA_DEFAULT = 1.0;
const double SIGMA_DEFAULT = 1.0;
const double REGULARIZATION_DEFAULT = 0;
const int TRAIN_MAX_ITER_DEFAULT = 1000;
const double HYPERPARAM_DATA_DEFAULT = 1;
const double KD_HYPERPARAM_DATA_DEFAULT = 1;
const double LIG_HYPERPARAM_DATA_DEFAULT = 1;
const double LIGAND_BONUS_DEFAULT = 90; // FMN

// function prototypes
void Usage(const Options &options);
void Version();
void ParseArguments(int argc, char **argv, Options &options, std::vector<std::string> &filenames);
void MakeFileDescriptions(const Options &options, const std::vector<std::string> &filenames, std::vector<FileDescription> &descriptions);

void RunGradientSanityCheck(const Options &options, const std::vector<FileDescription> &descriptions);

void RunEnergyTest(const Options &options, const std::vector<FileDescription> &descriptions);

void RunTrainingMode(const Options &options, const std::vector<FileDescription> &descriptions);

void RunPredictionMode(const Options &options, const std::vector<FileDescription> &descriptions);

void RunSampleMode(const Options &options, const std::vector<FileDescription> &descriptions);

void RunREVIMode(const Options &options, const std::vector<FileDescription> &descriptions);


void RunPredictionFoldChangeMode(const Options &options, const std::vector<FileDescription> &descriptions);

// default parameters
#include "Defaults.hpp"

/////////////////////////////////////////////////////////////////
// main()
//
// Main program.
/////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
#ifdef MULTI
//std::cout << "Initializing MPI..."<< std::endl;
    MPI_Init(&argc, &argv);
#endif

    // first, parse arguments
    Options options;
    std::vector<std::string> filenames;
    ParseArguments(argc, argv, options, filenames);

    // print out arguments to see that we parsed them correctly
    std::cerr << "Training mode: " << options.GetStringValue("training_mode") << std::endl
              << "Use constraints: " << options.GetBoolValue("use_constraints") << std::endl
              << "Use evidence: " << options.GetBoolValue("use_evidence") << std::endl;

    // second, read input files
    std::vector<FileDescription> descriptions;
    MakeFileDescriptions(options, filenames, descriptions);

/// check files
//     for (size_t i = 0; i < filenames.size(); i++)
//     {
//             const SStruct &sstruct = descriptions[i].sstruct.GetKd;
//             std::cout << "checking kd" << sstruct.GetKdData() << std::endl;
//     }
///

    // perform required task
    if (options.GetBoolValue("gradient_sanity_check"))
    {
        RunGradientSanityCheck(options, descriptions);
    }
    else if (options.GetBoolValue("test_energies"))
    {
        RunEnergyTest(options, descriptions);
    }

    else if (options.GetStringValue("training_mode") != "")
    {
        RunTrainingMode(options, descriptions);
    }
    else if (options.GetBoolValue("run_sample_mode")){
      RunSampleMode(options, descriptions);
    }
    else if (options.GetBoolValue("run_revi_mode")){
      RunREVIMode(options, descriptions);
    }
    else if (options.GetStringValue("prediction_mode") == "foldchange")
    {
        RunPredictionFoldChangeMode(options, descriptions);
    }
    else
    {
        RunPredictionMode(options, descriptions);
    }
    
#ifdef MULTI
    MPI_Finalize();
#endif
}

/////////////////////////////////////////////////////////////////
// Usage()
//
// Display program usage.
/////////////////////////////////////////////////////////////////

void Usage(const Options &options)
{
    std::cerr << std::endl
              << "Usage: contrafold [predict|predict-foldchange|sample|revi|train|em-train] [OPTION]... INFILE(s)" << std::endl 
              << std::endl
              << "       where [OPTION]...   is a list of zero or more optional arguments" << std::endl
              << "             INFILE(s)     is the name of the input BPSEQ, plain text, or FASTA file(s)" << std::endl 
              << std::endl
              << "Miscellaneous arguments:" << std::endl
              << "  --version                display program version information" << std::endl
              << "  --verbose                show detailed console output" << std::endl
              << "  --logbase LOG_BASE       set base of log-sum-exp" << std::endl
              << "  --viterbi                use Viterbi instead of posterior decoding for prediction, " << std::endl
              << "                           or max-margin instead of log-likelihood for training" << std::endl
              << "  --noncomplementary       allow non-{AU,CG,GU} pairs" << std::endl
              << std::endl 
              << "Additional arguments for 'predict' mode:" << std::endl
              << "  --params FILENAME        use particular model parameters" << std::endl
              << "  --constraints            use existing constraints (requires BPSEQ or FASTA format input)" << std::endl
              << "  --evidence               use experimental evidence (requires BPSEQ format input)" << std::endl
              << "  --centroid               use centroid estimator (as opposed to MEA estimator)" << std::endl
              << "  --gamma GAMMA            set sensivity/specificity tradeoff parameter (default: GAMMA=" << options.GetRealValue("gamma") << ")" << std::endl
              << "                             if GAMMA > 1, emphasize sensitivity" << std::endl
              << "                             if 0 <= GAMMA <= 1, emphasize specificity" << std::endl
              << "                             if GAMMA < 0, try tradeoff parameters of 2^k for k = -5,...,10" << std::endl
              << std::endl
              << "  --parens OUTFILEORDIR    write parenthesized output to file or directory" << std::endl
              << "  --bpseq OUTFILEORDIR     write BPSEQ output to file or directory" << std::endl
              << "  --posteriors CUTOFF OUTFILEORDIR" << std::endl
              << "                           write posterior pairing probabilities to file or directory" << std::endl
              << "  --partition              compute the partition function or Viterbi score only" << std::endl
              << std::endl
              << "Additional arguments for 'sample' mode:" <<std::endl
              << "  --nsamples N             number of samples (default 100)" << std::endl
              << "  --kappa k                weight to place on chemical mapping data (default 1.0)" << std::endl
              << std::endl
              << "Additional arguments for 'revi' mode:" <<std::endl
              << "  --sigma S                weight to place on data vs. energy model" << std::endl
              << std::endl
              << "Additional arguments for training (many input files may be specified):" << std::endl
              << "  --examplefile            read list of input files from provided text file (instead of as arguments)" << std::endl
              << "  --ligand                 include ligand MS2 kd prediction for riboswitches [HKWS]" << std::endl
              << "  --sanity                 perform gradient sanity check" << std::endl
              << "  --holdout F              use fraction F of training data for holdout cross-validation" << std::endl
              << "  --init_reg_file          use existing regularization weights, in file, to continue holdout CV OPT2 training" << std::endl
              << "  --regularize C           perform BFGS training, using a single regularization coefficient C" << std::endl
              << "  --maxiter N              for single regularization coefficient the max number of iterations" << std::endl
              << "  --shuffle                If included, shuffle data from datalist" << std::endl

              << "  --hyperparam_data K      weight on [chemical mapping] data-only examples" << std::endl
              << "  --initweights w          for single regularization coefficient an initial set of weights" << std::endl
              << "  --numdatasources n       the number of data sources for em-train" << std::endl
              << "  --batchsize b            mini-batch size for stochastic gradient training" << std::endl
              << "  --s0 s0                  stepsize for SGD is s0/(1+iter)^s1" << std::endl
              << "  --s1 s1                  stepsize for SGD is s0/(1+iter)^s1" << std::endl
              << std::endl
              << "Arguments added for EternaFold riboswitch training:" << std::endl

              << "  --kd_hyperparam_data K   weight on kd data examples" << std::endl
              << "  --lig_hyperparam_data K  weight on kd ligand data examples" <<std::endl
              << "  --ligand_bonus K         ligand bonus used in ligand-bound MS2 kd estimation" <<std::endl

              << std::endl;
    exit(0);
}

/////////////////////////////////////////////////////////////////
// Version()
//
// Display program version.
/////////////////////////////////////////////////////////////////

void Version()
{
#if PROFILE
    std::cerr << "CONTRAfold-SE(m) version 1 (NOT CURRENTLY IN USE) - Multiple sequence RNA secondary structure prediction based on CONTRAfold version 2.02" << std::endl << std::endl
#else
    std::cerr << "CONTRAfold-SE version 1 - RNA secondary structure prediction based on CONTRAfold version 2.02" << std::endl << std::endl
#endif
              << "Written by Cristina Pop and Chuan-Sheng Foo based on code by Chuong B. Do" << std::endl;
    exit(0);
}

/////////////////////////////////////////////////////////////////
// ParseArguments()
//
// Parse command line parameters.
/////////////////////////////////////////////////////////////////

void ParseArguments(int argc,
                    char **argv,
                    Options &options,
                    std::vector<std::string> &filenames)
{
    // register default options
    options.SetStringValue("training_mode", "");
    options.SetBoolValue("run_sample_mode",false);
    options.SetBoolValue("run_revi_mode",false);
    options.SetStringValue("prediction_mode", "");

    options.SetBoolValue("verbose_output", false);
    options.SetRealValue("log_base", 1.0);
    options.SetBoolValue("viterbi_parsing", false);
    options.SetBoolValue("train_with_ligand_data", false);
    options.SetBoolValue("allow_noncomplementary", false);

    options.SetStringValue("parameter_filename", "");
    options.SetBoolValue("use_constraints", false);
    options.SetBoolValue("centroid_estimator", false);
    options.SetBoolValue("use_evidence", false);
    options.SetRealValue("gamma", GAMMA_DEFAULT);
    options.SetRealValue("kappa", KAPPA_DEFAULT);
    options.SetRealValue("sigma", SIGMA_DEFAULT);
    options.SetStringValue("output_parens_destination", "");
    options.SetStringValue("output_bpseq_destination", "");
    options.SetRealValue("output_posteriors_cutoff", 0);
    options.SetStringValue("output_posteriors_destination", "");
    options.SetBoolValue("partition_function_only", false);
    options.SetBoolValue("shuffle", false);

    options.SetRealValue("HKWS_gradient_mode", 0);
    options.SetBoolValue("gradient_sanity_check", false);
    options.SetBoolValue("test_energies", false);
    options.SetRealValue("holdout_ratio", 0);
    options.SetStringValue("opt2_regularization_weights", "");

    options.SetRealValue("regularization_coefficient", REGULARIZATION_DEFAULT);
    options.SetRealValue("ligand_bonus", LIGAND_BONUS_DEFAULT);

    options.SetStringValue("train_examplefile", "");
    options.SetIntValue("train_max_iter", TRAIN_MAX_ITER_DEFAULT);
    options.SetStringValue("train_initweights_filename", "");
    options.SetStringValue("train_priorweights_filename", "");
    options.SetIntValue("num_data_sources",0);
    options.SetIntValue("batch_size", 1);
    options.SetRealValue("s0", 0.0001);
    options.SetRealValue("s1", 0);
    options.SetRealValue("hyperparam_data",HYPERPARAM_DATA_DEFAULT);
    options.SetRealValue("kd_hyperparam_data",KD_HYPERPARAM_DATA_DEFAULT);
    options.SetRealValue("lig_hyperparam_data",LIG_HYPERPARAM_DATA_DEFAULT);
    options.SetIntValue("nsamples",100);


    // check for sufficient arguments
    if (argc < 2) Usage(options);
    filenames.clear();

    // check for prediction or training mode    
    if (!strcmp(argv[1], "train")) 
    {
        options.SetStringValue("training_mode", "supervised");
    }
    else if (!strcmp(argv[1], "em-train"))
    {
        options.SetStringValue("training_mode", "em");
    }
    else if (!strcmp(argv[1], "em-train-sgd"))
    {
        options.SetStringValue("training_mode", "em-sgd");
    }
    else if (!strcmp(argv[1], "sample"))
    {
        options.SetBoolValue("run_sample_mode",true);
    }
    else if (!strcmp(argv[1], "revi"))
    {
        options.SetBoolValue("run_revi_mode",true);
    }
    else if (!strcmp(argv[1], "predict-foldchange"))
    {
        options.SetStringValue("prediction_mode", "foldchange");
    }
    else if (strcmp(argv[1], "predict"))
    {
        Error("CONTRAfold must be run in either 'predict', 'train', 'sample',or 'revi' mode.");
    }

    // go through remaining arguments
    for (int argno = 2; argno < argc; argno++)
    {
        // parse optional arguments
        if (argv[argno][0] == '-')
        {
            // miscellaneous options
            if (!strcmp(argv[argno], "--version"))
            {
                Version();
            }
            else if (!strcmp(argv[argno], "--verbose"))
            {
                options.SetBoolValue("verbose_output", true);
            }
            else if (!strcmp(argv[argno], "--logbase"))
            {
                if (argno == argc - 1) Error("Must specify log base LOG_BASE after --logbase.");
                double value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse log base.");
                if (value <= 0)
                    Error("Log base must be positive.");
                options.SetRealValue("log_base", value);
            }
            else if (!strcmp(argv[argno], "--viterbi"))
            {
                options.SetBoolValue("viterbi_parsing", true);
            }
            else if (!strcmp(argv[argno], "--noncomplementary"))
            {
                options.SetBoolValue("allow_noncomplementary", true);
            }
            
            // prediction options
            else if (!strcmp(argv[argno], "--params"))
            {
                if (argno == argc - 1) Error("Must specify FILENAME after --params.");
                options.SetStringValue("parameter_filename", argv[++argno]);
            }
            else if (!strcmp(argv[argno], "--constraints"))
            {
                options.SetBoolValue("use_constraints", true);
            }
            else if (!strcmp(argv[argno], "--evidence"))
            {
                options.SetBoolValue("use_evidence", true);
            }
            else if (!strcmp(argv[argno], "--centroid"))
            {
                options.SetBoolValue("centroid_estimator", true);
            }
            else if (!strcmp(argv[argno], "--gamma"))
            {
                if (argno == argc - 1) Error("Must specify trade-off parameter GAMMA after --gamma.");
                double value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse value after --gamma.");
                options.SetRealValue("gamma", value);
            }
            else if (!strcmp(argv[argno], "--kappa"))
            {
                if (argno == argc - 1) Error("Must specify chemical mapping reweighting parameter KAPPA after --kappa.");
                double value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse value after --kappa.");
                options.SetRealValue("kappa", value);
            }
            else if (!strcmp(argv[argno], "--sigma"))
            {
                if (argno == argc - 1) Error("Must specify data reweighting parameter SIGMA after --sigma.");
                double value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse value after --sigma.");
                options.SetRealValue("sigma", value);
            }
            else if (!strcmp(argv[argno], "--parens"))
            {
                if (argno == argc - 1) Error("Must specify output file or directory name after --parens.");
                options.SetStringValue("output_parens_destination", argv[++argno]);
            }
            else if (!strcmp(argv[argno], "--bpseq"))
            {
                if (argno == argc - 1) Error("Must specify output file or directory name after --bpseq.");
                options.SetStringValue("output_bpseq_destination", argv[++argno]);
            }
            else if (!strcmp(argv[argno], "--posteriors"))
            {
                if (argno == argc - 1) Error("Must specify posterior probability threshold CUTOFF after --posteriors.");
                double value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse cutoff value after --posteriors.");
                options.SetRealValue("output_posteriors_cutoff", value);
                if (argno == argc - 1) Error("Must specify output file or directory for --posteriors.");
                options.SetStringValue("output_posteriors_destination", argv[++argno]);
            }
            else if (!strcmp(argv[argno], "--partition"))
            {
                options.SetBoolValue("partition_function_only", true);
            }
            else if (!strcmp(argv[argno], "--shuffle"))
            {
                options.SetBoolValue("shuffle", true);
            }
            // training options
            else if (!strcmp(argv[argno], "--examplefile"))
            {
                if (argno == argc - 1) Error("Must specify filename after --examplefile.");
                options.SetStringValue("train_examplefile", argv[++argno]);
            }
            else if (!strcmp(argv[argno], "--sanity"))
            {
                options.SetBoolValue("gradient_sanity_check", true);
            }
            else if (!strcmp(argv[argno], "--test_energies"))
            {
                options.SetBoolValue("test_energies", true);
            }
            else if (!strcmp(argv[argno], "--ligand"))
            {
                options.SetBoolValue("train_with_ligand_data", true);
            }
            else if (!strcmp(argv[argno], "--ligand_bonus"))
            {
                if (argno == argc - 1) Error("Must specify ligand bonus value after --ligand_bonus.");
                double value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse regularization parameter after --ligand_bonus.");
                if (value < 0)
                    Error("ligand bonus should not be negative.");
                options.SetRealValue("ligand_bonus", value);
            }
            else if (!strcmp(argv[argno], "--HKWS_gradient_mode"))
            {
                if (argno == argc - 1) Error("Must specify gradient mode.");
                double value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse --HKWS_gradient_mode.");
                if (value < 0)
                    Error("HKWS_gradient_mode should not be negative.");
                options.SetRealValue("HKWS_gradient_mode", value);  
            }
            else if (!strcmp(argv[argno], "--holdout"))
            {
                if (argno == argc - 1) Error("Must specify holdout ratio F after --holdout.");
                double value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse holdout ratio.");
                if (value < 0 || value > 1)
                    Error("Holdout ratio must be between 0 and 1.");
                options.SetRealValue("holdout_ratio", value);
            }
            else if (!strcmp(argv[argno], "--regularize"))
            {
                if (argno == argc - 1) Error("Must specify regularization parameter C after --regularize.");
                double value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse regularization parameter after --regularize.");
                if (value < 0)
                    Error("Regularization parameter should not be negative.");
                options.SetRealValue("regularization_coefficient", value);
            }

            else if (!strcmp(argv[argno], "--init_reg_file"))
            {
                if (argno == argc - 1) Error("Must specify filename containing OPT2 weights after --init_reg_file.");
                options.SetStringValue("opt2_regularization_weights", argv[++argno]);
            }

            else if (!strcmp(argv[argno], "--maxiter"))
            {
                if (argno == argc - 1) Error("Must specify max number of iterations after --maxiter.");
                int value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse max number of iterations after --maxiter.");
                if (value < 0)
                    Error("Max number of iterations should not be negative.");
                options.SetIntValue("train_max_iter", value);
            }

            else if (!strcmp(argv[argno], "--nsamples"))
            {
                if (argno == argc - 1) Error("Must specify number of samples after --nsamples.");
                int value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse max number of iterations after --nsamples.");
                if (value < 0)
                    Error("Number of samples should not be negative.");
                options.SetIntValue("nsamples", value);
            }
            else if (!strcmp(argv[argno], "--hyperparam_data"))
            {
                if (argno == argc - 1) Error("Must specify a value after --hyperparam_data.");
                double value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse value after --hyperparam_data.");
                if (value < 0)
                    Error("Value after --hyperparam_data should not be negative.");
                options.SetRealValue("hyperparam_data", value);
            }
            else if (!strcmp(argv[argno], "--kd_hyperparam_data"))
            {
                if (argno == argc - 1) Error("Must specify a value after --kd_hyperparam_data.");
                double value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse value after --kd_hyperparam_data.");
                if (value < 0)
                    Error("Value after --kd_hyperparam_data should not be negative.");
                options.SetRealValue("kd_hyperparam_data", value);
            }
            else if (!strcmp(argv[argno], "--lig_hyperparam_data"))
            {
                if (argno == argc - 1) Error("Must specify a value after --lig_hyperparam_data.");
                double value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse value after --lig_hyperparam_data.");
                if (value < 0)
                    Error("Value after --lig_hyperparam_data should not be negative.");
                options.SetRealValue("lig_hyperparam_data", value);
            }
            else if (!strcmp(argv[argno], "--initweights"))
            {
                if (argno == argc - 1) Error("Must specify filename of initial weights after --initweights.");
                options.SetStringValue("train_initweights_filename", argv[++argno]);
            }
            else if (!strcmp(argv[argno], "--priorweights"))
            {
                if (argno == argc - 1) Error("Must specify filename of prior weights after --priorweights.");
                options.SetStringValue("train_priorweights_filename", argv[++argno]);
            }
            else if (!strcmp(argv[argno],"--numdatasources"))
            {
                if (argno == argc - 1) Error("Must specific an integer number of data sources after --numdatasources.");
                int value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse number of data sources after --numdatasources.");
                if (value < 0)
                    Error("Number of data sources should be positive: %d", value);
                if (value == 0)
                    std::cout << "Warning: number of data sources set to zero!" << std::endl;
                options.SetIntValue("num_data_sources",value);
            }
            else if (!strcmp(argv[argno], "--batchsize"))
            {
                if (argno == argc - 1) Error("Must specify an integer batch size after --batchsize.");
                int value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse number of data sources after --batchsize.");
                if (value <= 0)
                    Error("Batch size should be positive: %d", value);
                options.SetIntValue("batch_size", value);
            }
            else if (!strcmp(argv[argno], "--s0"))
            {
                if (argno == argc - 1) Error("Must specify stepsize parameter s0 after --s0");
                double value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse stepsize parameter after --s0.");
                if (value < 0)
                    Error("Stepsize parameter should not be negative.");
                options.SetRealValue("s0", value);
            }
            else if (!strcmp(argv[argno], "--s1"))
            {
                if (argno == argc - 1) Error("Must specify stepsize parameter s1 after --s1");
                double value;
                if (!ConvertToNumber(argv[++argno], value))
                    Error("Unable to parse stepsize parameter after --s1.");
                if (value < 0)
                    Error("Stepsize parameter should not be negative.");
                options.SetRealValue("s1", value);
            }
            else
            {
                Error("Unknown option \"%s\" specified.  Run program without any arguments to see command-line options.", argv[argno]);
            }
        }
        else
        {
            filenames.push_back(argv[argno]);
        }
    }

    // ensure that at least one input file specified
    if (filenames.size() == 0 && options.GetStringValue("train_examplefile") == "")
        Error("No filenames / example file specified.");

    if (filenames.size() != 0 && options.GetStringValue("train_examplefile") != "")
        Error("Filenames specified on command line and via --examplefile. Use either.");

    if (options.GetStringValue("train_examplefile") != "")
    {
        std::ifstream examplefile(options.GetStringValue("train_examplefile").c_str());
        if (examplefile.fail()) Error("Unable to open training example file: %s", options.GetStringValue("train_examplefile").c_str());

        std::string s;
        while (std::getline(examplefile, s))
        {
            s = Trim(s);
            filenames.push_back(s);
        }

        if (filenames.size() == 0)
            Error("No files read from train example file!");

        std::cerr << "Read " << filenames.size() << " files from " << options.GetStringValue("train_examplefile") << std::endl;
    }

    // check to make sure that arguments make sense
    if (options.GetStringValue("training_mode") == "")
    {
        if (options.GetStringValue("train_initweights_filename") != "")
            Error("The --initweights flag is not used outside of train mode.");
        if (options.GetRealValue("hyperparam_data") != HYPERPARAM_DATA_DEFAULT)
            Error("The --hyperparam_data flag is not used outside of training mode.");
        if (options.GetRealValue("kd_hyperparam_data") != KD_HYPERPARAM_DATA_DEFAULT)
            Error("The --kd_hyperparam_data flag is not used outside of training mode.");
        if (options.GetRealValue("lig_hyperparam_data") != LIG_HYPERPARAM_DATA_DEFAULT)
            Error("The --lig_hyperparam_data flag is not used outside of training mode.");
    }

    // check to make sure that arguments make sense
    if (options.GetStringValue("training_mode") != "em" &&
        options.GetStringValue("training_mode") != "em-sgd" &&
        options.GetStringValue("training_mode") != "supervised")
    {
        if (options.GetIntValue("train_max_iter") != TRAIN_MAX_ITER_DEFAULT)
            Error("The --maxiter flag is not used outside of training mode.");
    }

    // check to make sure that arguments make sense
    if (options.GetStringValue("training_mode") == "em")
    {
        if (options.GetIntValue("num_data_sources") <= 0)
            Error("Registered <= 0 number of data sources. Please verify this is correct.");
    }
    
    // check to make sure that arguments make sense
    if (options.GetStringValue("training_mode") != "")
    {
        if (options.GetStringValue("parameter_filename") != "")
            Error("Should not specify parameter file for training mode.");
        if (options.GetBoolValue("use_constraints") && options.GetStringValue("training_mode") != "em")
            Error("The --constraints flag has no effect in training mode.");
        if (options.GetBoolValue("centroid_estimator"))
            Error("The --centroid flag has no effect in training mode.");
        if (options.GetRealValue("gamma") != GAMMA_DEFAULT)
            Error("Gamma parameter should not be specified in training mode.");
        if (options.GetStringValue("output_parens_destination") != "")
            Error("The --parens option cannot be used in training mode.");
        if (options.GetStringValue("output_bpseq_destination") != "")
            Error("The --bpseq option cannot be used in training mode.");
        if (options.GetStringValue("output_posteriors_destination") != "" ||
            options.GetRealValue("output_posteriors_cutoff") != 0)
            Error("The --posteriors option cannot be used in training mode.");
        if (options.GetBoolValue("partition_function_only"))
            Error("The --partition flag cannot be used in training mode.");
        if (options.GetRealValue("regularization_coefficient") != REGULARIZATION_DEFAULT &&
            options.GetRealValue("holdout_ratio") > 0)
            Error("The --holdout and --regularize options cannot be specified simultaneously.");
    }
    else
    {
        if (options.GetRealValue("gamma") < 0 &&
            options.GetStringValue("output_parens_destination") == "" &&
            options.GetStringValue("output_bpseq_destination") == "" &&
            options.GetStringValue("output_posteriors_destination") == "")
            Error("Output directory must be specified when using GAMMA < 0.");
    
        if (options.GetBoolValue("use_constraints") && options.GetBoolValue("use_evidence"))
            Error("You can only use either constraints or evidence, not both together.");

#ifdef MULTI
        if (filenames.size() > 1 &&
            options.GetStringValue("output_parens_destination") == "" &&
            options.GetStringValue("output_bpseq_destination") == "" &&
            options.GetStringValue("output_posteriors_destination") == "")
            Error("Output directory must be specified when performing predictions for multiple input files.");
#endif
        if (options.GetBoolValue("viterbi_parsing") &&
            options.GetStringValue("output_posteriors_destination") != "")
            Error("The --posteriors option cannot be used with Viterbi parsing.");
    }
}

/////////////////////////////////////////////////////////////////
// MakeFileDescriptions()
//
// Build file descriptions
/////////////////////////////////////////////////////////////////

void MakeFileDescriptions(const Options &options,
                          const std::vector<std::string> &filenames,
                          std::vector<FileDescription> &descriptions)
{
    descriptions.clear();
    for (size_t i = 0; i < filenames.size(); i++)
    {
        descriptions.push_back(FileDescription(filenames[i],
                                               options.GetBoolValue("allow_noncomplementary"),options.GetIntValue("num_data_sources")));
    }

    std::sort(descriptions.begin(), descriptions.end());

}
/////////////////////////////////////////////////////////////////
// RunEnergyTest()
//
// Print energies for test sequence.
/////////////////////////////////////////////////////////////////

void RunEnergyTest(const Options &options,
                            const std::vector<FileDescription> &descriptions)
{
    ParameterManager parameter_manager;
    InferenceEngine inference_engine(options.GetBoolValue("allow_noncomplementary"),
        options.GetIntValue("num_data_sources"), options.GetRealValue("kappa"));
    inference_engine.RegisterParameters(parameter_manager);
    ComputationEngine computation_engine(options, descriptions, inference_engine, parameter_manager);
    ComputationWrapper computation_wrapper(computation_engine);
    
    // decide whether I'm a compute node or master node
    if (computation_engine.IsComputeNode())
    {
        computation_engine.RunAsComputeNode();
        return;
    }

    const std::string output_parens_destination = options.GetStringValue("output_parens_destination");
    const std::string output_bpseq_destination = options.GetStringValue("output_bpseq_destination");
    const std::string output_posteriors_destination = options.GetStringValue("output_posteriors_destination");

    // load parameters
    std::vector<RealT> w;

    if (options.GetStringValue("parameter_filename") != "")
    {
        parameter_manager.ReadFromFile(options.GetStringValue("parameter_filename"), w);
    }
    else
    {
#if PROFILE
        w = GetDefaultProfileValues();
#else
        if (options.GetBoolValue("allow_noncomplementary"))
            w = GetDefaultNoncomplementaryValues();
        else
            w = GetDefaultComplementaryValues();
#endif
    }

    
        // create output directories for output files, if needed

    computation_wrapper.TestEnergies(computation_wrapper.GetAllUnits(), w, options.GetRealValue("gamma"), options.GetRealValue("log_base"));
    
    computation_engine.StopComputeNodes();
}
/////////////////////////////////////////////////////////////////
// RunGradientSanityCheck()
//
// Compute gradient sanity check.
/////////////////////////////////////////////////////////////////

void RunGradientSanityCheck(const Options &options,
                            const std::vector<FileDescription> &descriptions)
{
    // The architecture of the code is somewhat complicated here, so
    // here's a quick explanation:
    // 
    //    ParameterManager: associates each parameter of the model
    //                      with a name and manages hyperparameter
    //                      groups
    //                     
    //    InferenceEngine: performs application-specific
    //                     (loss-augmented) inference
    //
    //    ComputationEngine: makes all necessary calls to dynamic
    //                       programming routines for processing
    //                       individual sequences and interfaces with
    //                       distributed computation module
    //
    //    ComputationWrapper: provides a high-level interface for
    //                        performing computations on groups of
    //                        sequences
    //
    //    OuterOptimizationWrapper / InnerOptimizationWrapper:
    //                        interface between computation routines
    //                        and optimization routines
    
    ParameterManager parameter_manager;
    InferenceEngine inference_engine(options.GetBoolValue("allow_noncomplementary"),
      options.GetIntValue("num_data_sources"), options.GetRealValue("kappa"));
    inference_engine.RegisterParameters(parameter_manager);
    ComputationEngine computation_engine(options, descriptions, inference_engine, parameter_manager);
    ComputationWrapper computation_wrapper(computation_engine);

    // decide whether I'm a compute node or master node
    if (computation_engine.IsComputeNode())
    {
        computation_engine.RunAsComputeNode();
        return;
    }

    std::vector<RealT> w;  // (parameter_manager.GetNumLogicalParameters(), RealT(0));
    const std::string initweights_filename = options.GetStringValue("train_initweights_filename");
    if (initweights_filename != "")
    {
        // if specified
        parameter_manager.ReadFromFile(initweights_filename, w);
    }
    else
    {
        // to 0 otherwise
        for (int i = 0; i < (int)parameter_manager.GetNumLogicalParameters(); i++)
            w.push_back(RealT(0));
    }


    computation_wrapper.SanityCheckGradient(computation_wrapper.GetAllUnits(), w);
    computation_engine.StopComputeNodes();
}

/////////////////////////////////////////////////////////////////
// RunTrainingMode()
//
// Run CONTRAfold in training mode.
/////////////////////////////////////////////////////////////////

void RunTrainingMode(const Options &options,
                     const std::vector<FileDescription> &descriptions)
{
    ParameterManager parameter_manager;
    InferenceEngine inference_engine(options.GetBoolValue("allow_noncomplementary"),
      options.GetIntValue("num_data_sources"), options.GetRealValue("kappa"));
    inference_engine.RegisterParameters(parameter_manager);
    ComputationEngine computation_engine(options, descriptions, inference_engine, parameter_manager);
    ComputationWrapper computation_wrapper(computation_engine);


    // decide whether I'm a compute node or master node
    if (computation_engine.IsComputeNode())
    {
        computation_engine.RunAsComputeNode();
        return;
    }

    const std::string initweights_filename = options.GetStringValue("train_initweights_filename");
    const std::string priorweights_filename = options.GetStringValue("train_priorweights_filename");
    const int train_max_iter = options.GetIntValue("train_max_iter");
    const RealT hyperparam_data = options.GetRealValue("hyperparam_data");
    const RealT kd_hyperparam_data = options.GetRealValue("kd_hyperparam_data");
    const RealT lig_hyperparam_data = options.GetRealValue("lig_hyperparam_data");


    // set the initial parameters:
    std::vector<RealT> w;
    if (initweights_filename != "")
    {
        // if specified
        parameter_manager.ReadFromFile(initweights_filename, w);
    }
    else
    {
        // to 0 otherwise
        for (int i = 0; i < (int)parameter_manager.GetNumLogicalParameters(); i++)
            w.push_back(RealT(0));
    }

    std::vector<RealT> w0;
    if (priorweights_filename != "")
    {
        // if specified
        parameter_manager.ReadFromFile(priorweights_filename, w0);
    }
    else
    {
        // to 0 otherwise
        for (int i = 0; i < (int)parameter_manager.GetNumLogicalParameters(); i++)
            w0.push_back(RealT(0));
    }

    std::vector<int> units = computation_wrapper.FilterNonparsable(computation_wrapper.GetAllUnits());

    OptimizationWrapper optimization_wrapper(computation_wrapper);
    
    std::vector<RealT> regularization_coefficients;

    // decide between using a fixed regularization parameter or
    // using cross-validation to determine regularization parameters
    if (options.GetRealValue("holdout_ratio") <= 0) {
        if (options.GetStringValue("opt2_regularization_weights") != ""){
            // if provided file of regularization weights, read in
            std::ifstream regfile(options.GetStringValue("opt2_regularization_weights").c_str());
            if (regfile.fail()) Error("Could not open regfile for reading.");

            std::string line;

            while (std::getline(regfile, line))
            {
                float value;
                std::stringstream ss(line);

                while (ss >> value)
                {
                    regularization_coefficients.push_back(value);
                }
            }

            //TODO: write check that length log_C = NumParameterGroups

        } else {
            std::cout << "no opt2_regularization_weights file, setting all to be same" << std::endl;
            //regularization_coefficients(parameter_manager.GetNumParameterGroups(), options.GetRealValue("regularization_coefficient"));

            for (size_t i = 0; i < parameter_manager.GetNumParameterGroups(); i++)
                regularization_coefficients.push_back(options.GetRealValue("regularization_coefficient"));
        }
        std::cout << "reg coefficients" << regularization_coefficients << std::endl;
        if (options.GetStringValue("training_mode") == "em") {
            optimization_wrapper.TrainEM(units, w, regularization_coefficients,train_max_iter);
        } else if (options.GetStringValue("training_mode") == "em-sgd") {
        // Don't regularize evidence CPD parameters
            std::cout << parameter_manager.GetNumParameterGroups() << std::endl;
            if (parameter_manager.GetNumParameterGroups() != 2)
                Error("Using em-sgd with multiple hyperparameters is not supported");
            regularization_coefficients[1] = 0;
            optimization_wrapper.TrainSGD(units, w, regularization_coefficients);

        } else {
                optimization_wrapper.Train(units, w, w0, regularization_coefficients);

        }
    }
    else
    {
        // TODO: add hyperparam_data for sgd mode
        // if (options.GetStringValue("training_mode") == "em") // HW commenting out because not defined in case of not hyperparam grid search
        //     optimization_wrapper.LearnHyperparametersEM(units, w, train_max_iter);
        // else
            optimization_wrapper.LearnHyperparameters(units, w);

    }
    
    parameter_manager.WriteToFile("optimize.params.final", w);
    computation_engine.StopComputeNodes();
}

/////////////////////////////////////////////////////////////////
// RunPredictionMode()
//
// Run CONTRAfold in prediction mode.
/////////////////////////////////////////////////////////////////

void RunPredictionMode(const Options &options,
                       const std::vector<FileDescription> &descriptions)
{
    ParameterManager parameter_manager;
    InferenceEngine inference_engine(options.GetBoolValue("allow_noncomplementary"),
      options.GetIntValue("num_data_sources"), options.GetRealValue("kappa"));
    inference_engine.RegisterParameters(parameter_manager);
    ComputationEngine computation_engine(options, descriptions, inference_engine, parameter_manager);
    ComputationWrapper computation_wrapper(computation_engine);
    
    // decide whether I'm a compute node or master node
    if (computation_engine.IsComputeNode())
    {
        computation_engine.RunAsComputeNode();
        return;
    }

    const std::string output_parens_destination = options.GetStringValue("output_parens_destination");
    const std::string output_bpseq_destination = options.GetStringValue("output_bpseq_destination");
    const std::string output_posteriors_destination = options.GetStringValue("output_posteriors_destination");

    // load parameters
    std::vector<RealT> w;

    if (options.GetStringValue("parameter_filename") != "")
    {
        parameter_manager.ReadFromFile(options.GetStringValue("parameter_filename"), w);
    }
    else
    {
#if PROFILE
        w = GetDefaultProfileValues();
#else
        if (options.GetBoolValue("allow_noncomplementary"))
            w = GetDefaultNoncomplementaryValues();
        else
            w = GetDefaultComplementaryValues();
#endif
    }

    if (options.GetRealValue("gamma") < 0)
    {
        // create directories for storing each run
        if (output_parens_destination != "") MakeDirectory(output_parens_destination);
        if (output_bpseq_destination != "") MakeDirectory(output_bpseq_destination);
        if (output_posteriors_destination != "") MakeDirectory(output_posteriors_destination);
        
        // try different values of gamma
        for (int k = -5; k <= 10; k++)
        {
            // create output subdirectories, if needed
            const double gamma = Pow(2.0, double(k));

            if (descriptions.size() > 1)
            {
                if (output_parens_destination != "")
                    MakeDirectory(SPrintF("%s%c%s.gamma=%lf",
                                          output_parens_destination.c_str(),
                                          DIR_SEPARATOR_CHAR,
                                          GetBaseName(output_parens_destination).c_str(), gamma));
                if (output_bpseq_destination != "")
                    MakeDirectory(SPrintF("%s%c%s.gamma=%lf",
                                          output_bpseq_destination.c_str(),
                                          DIR_SEPARATOR_CHAR,
                                          GetBaseName(output_bpseq_destination).c_str(), gamma));
                if (output_posteriors_destination != "")
                    MakeDirectory(SPrintF("%s%c%s.gamma=%lf",
                                          output_posteriors_destination.c_str(),
                                          DIR_SEPARATOR_CHAR,
                                          GetBaseName(output_posteriors_destination).c_str(), gamma));
            }
            
            // perform predictions
            computation_wrapper.Predict(computation_wrapper.GetAllUnits(), w, gamma, options.GetRealValue("log_base"));
        }
    }
    else
    {
        // create output directories for output files, if needed
        if (descriptions.size() > 1)
        {
            if (output_parens_destination != "") MakeDirectory(output_parens_destination);
            if (output_bpseq_destination != "") MakeDirectory(output_bpseq_destination);
            if (output_posteriors_destination != "") MakeDirectory(output_posteriors_destination);
        }
        
        computation_wrapper.Predict(computation_wrapper.GetAllUnits(), w, options.GetRealValue("gamma"), options.GetRealValue("log_base"));
    }
    computation_engine.StopComputeNodes();
}

/////////////////////////////////////////////////////////////////
// RunPredictionFoldChangeMode() HKWS
//
// Run CONTRAfold in prediction foldchange mode.
/////////////////////////////////////////////////////////////////

void RunPredictionFoldChangeMode(const Options &options,
                       const std::vector<FileDescription> &descriptions)
{
    ParameterManager parameter_manager;
    InferenceEngine inference_engine(options.GetBoolValue("allow_noncomplementary"),
      options.GetIntValue("num_data_sources"), options.GetRealValue("kappa"));
    inference_engine.RegisterParameters(parameter_manager);
    ComputationEngine computation_engine(options, descriptions, inference_engine, parameter_manager);
    ComputationWrapper computation_wrapper(computation_engine);
    //std::cout << "unconditional_score,conditional_score,conditional_score2,conditional_score3,pred_log_kd_no_lig,true_log_kd_no_lig,pred_log_kd_w_lig,true_log_kd_no_lig" << std::endl;   
    // decide whether I'm a compute node or master node
    if (computation_engine.IsComputeNode())
    {
        computation_engine.RunAsComputeNode();
        return;
    }

    // const std::string output_parens_destination = options.GetStringValue("output_parens_destination");
    const std::string output_bpseq_destination = options.GetStringValue("output_bpseq_destination"); //TODO HKWS: make new var output_foldx_destination
    // const std::string output_posteriors_destination = options.GetStringValue("output_posteriors_destination");

   if (output_bpseq_destination != "")
   MakeDirectory(SPrintF("%s%c%s.gamma=%lf",
                         output_bpseq_destination.c_str(),
                         DIR_SEPARATOR_CHAR,
                         GetBaseName(output_bpseq_destination).c_str(), options.GetRealValue("gamma")));

    // load parameters
    std::vector<RealT> w;
    if (options.GetStringValue("parameter_filename") != "")
    {
        parameter_manager.ReadFromFile(options.GetStringValue("parameter_filename"), w);
    }
    else
    {
#if PROFILE
        w = GetDefaultProfileValues();
#else
        if (options.GetBoolValue("allow_noncomplementary"))
            w = GetDefaultNoncomplementaryValues();
        else
            w = GetDefaultComplementaryValues();
#endif
    }
    computation_wrapper.PredictFoldChange(computation_wrapper.GetAllUnits(), w, options.GetRealValue("gamma"), options.GetRealValue("log_base"));
    computation_engine.StopComputeNodes();
}

/////////////////////////////////////////////////////////////////
// RunSampleMode()
//
// Run CONTRAfold in sample mode.
/////////////////////////////////////////////////////////////////

void RunSampleMode(const Options &options,
                       const std::vector<FileDescription> &descriptions)
{
    ParameterManager parameter_manager;
    InferenceEngine inference_engine(options.GetBoolValue("allow_noncomplementary"),
      options.GetIntValue("num_data_sources"), options.GetRealValue("kappa"));
    inference_engine.RegisterParameters(parameter_manager);
    ComputationEngine computation_engine(options, descriptions, inference_engine, parameter_manager);
    ComputationWrapper computation_wrapper(computation_engine);
    
    // decide whether I'm a compute node or master node
    if (computation_engine.IsComputeNode())
    {
        computation_engine.RunAsComputeNode();
        return;
    }

    // const std::string output_parens_destination = options.GetStringValue("output_parens_destination");
    // const std::string output_bpseq_destination = options.GetStringValue("output_bpseq_destination");
    // const std::string output_posteriors_destination = options.GetStringValue("output_posteriors_destination");

    // load parameters
    std::vector<RealT> w;

    if (options.GetStringValue("parameter_filename") != "")
    {
        parameter_manager.ReadFromFile(options.GetStringValue("parameter_filename"), w);
    }
    else
    {
#if PROFILE
        w = GetDefaultProfileValues();
#else
        if (options.GetBoolValue("allow_noncomplementary"))
            w = GetDefaultNoncomplementaryValues();
        else
            w = GetDefaultComplementaryValues();
#endif
    }

        
    computation_wrapper.Sample(computation_wrapper.GetAllUnits(), w, options.GetRealValue("gamma"), options.GetRealValue("log_base"));
    
    computation_engine.StopComputeNodes();
}

/////////////////////////////////////////////////////////////////
// RunREVIMode()
//
// Run REVI.
/////////////////////////////////////////////////////////////////

void RunREVIMode(const Options &options,
                       const std::vector<FileDescription> &descriptions)
{
    ParameterManager parameter_manager;
    InferenceEngine inference_engine(options.GetBoolValue("allow_noncomplementary"),options.GetIntValue("num_data_sources"), options.GetRealValue("kappa"));
    inference_engine.RegisterParameters(parameter_manager);
    ComputationEngine computation_engine(options, descriptions, inference_engine, parameter_manager);
    ComputationWrapper computation_wrapper(computation_engine);
    
    // decide whether I'm a compute node or master node
    if (computation_engine.IsComputeNode())
    {
        computation_engine.RunAsComputeNode();
        return;
    }

    // load parameters
    std::vector<RealT> w;

    if (options.GetStringValue("parameter_filename") != "")
    {
        parameter_manager.ReadFromFile(options.GetStringValue("parameter_filename"), w);
    }
    else
    {
#if PROFILE
        w = GetDefaultProfileValues();
#else
        if (options.GetBoolValue("allow_noncomplementary"))
            w = GetDefaultNoncomplementaryValues();
        else
            w = GetDefaultComplementaryValues();
#endif
    }

    computation_wrapper.RunREVI(computation_wrapper.GetAllUnits(), w, options.GetRealValue("gamma"), options.GetRealValue("log_base"), options.GetRealValue("sigma"));
    
    computation_engine.StopComputeNodes();
}
