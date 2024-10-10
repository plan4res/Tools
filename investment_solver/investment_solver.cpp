/*--------------------------------------------------------------------------*/
/*----------------------- File investment_solver.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * This is a convenient tool for solving the investment problem defined by an
 * InvestmentBlock. The description of the InvestmentBlock must be given in a
 * netCDF file. This tool can be executed as follows:
 *
 *   ./investment_solver [-s] [-e] [-o] [-l FILE] [-n NUMBER] [-B FILE]
 *                       [-p PATH] [-c PATH] [-x FILE ] -S FILE <nc4-file>
 *
 * The only mandatory arguments are the netCDF file containing the description
 * of the InvestmentBlock and the solver configuration file indicated by the
 * -S option. This netCDF file can be either a BlockFile or a ProbFile. The
 * BlockFile can contain any number of child groups, each one describing an
 * InvestmentBlock. Every InvestmentBlock is then solved. The ProbFile can
 * also contain any number of child groups, each one having the description of
 * an InvestmentBlock alongside the description of a BlockConfig and a
 * BlockSolverConfig for the InvestmentBlock. Also in this case, every
 * InvestmentBlock is solved.
 *
 * The -c option specifies the prefix to the paths to all configuration
 * files. This means that if PATH is the value passed to the -c option, then
 * the name (or path) to each configuration file will be prepended by
 * PATH. The -p option specifies the prefix to the paths to all files
 * specified by the attribute "filename" in the input netCDF file.
 *
 * It is possible to provide an initial point (initial solution or initial
 * investment) through the -x option. This option must be followed by a file
 * containing the initial point. If there are N assets subject to investment,
 * then this file must contain N numbers, where the i-th number is the initial
 * value for the investment in the i-th asset. If this option is not used,
 * then the initial value x_i for the investment in the i-th asset is
 * determined as follows. If the lower bound l_i on the i-th investment is
 * finite, then x_i = l_i. Otherwise, if the upper bound u_i on the i-th
 * investment is finite, then x_i = u_i. Otherwise, if both bounds are not
 * finite, then x_i = 0.
 *
 * To simulate a given investment, i.e., to compute the investment function at
 * a given point, the -s option must be used. The investment to be simulated
 * is given by the initial point as described above: a given point provided by
 * the -x option or the default initial point.
 *
 * If the -o option is used, then part of the primal and dual solutions of
 * every UCBlock for each scenario is output while the investment function is
 * computed. Typically, one may want the solutions to be output in simulation
 * mode (i.e., when the -s option is used).
 *
 * The -n option specifies the number of sub-Blocks of SDDPBlock that must be
 * constructed for each stage. By default, SDDPBlock contains a single
 * sub-Blocks for each stage. This option must be provided in order to solve
 * multiple scenarios in parallel. In this case, the number of scenarios that
 * are solved in parallel is n (assuming n is not larger than the number of
 * scenarios).
 *
 * The -B and -S options are only considered if the given netCDF file is a
 * BlockFile. The -B option specifies a BlockConfig file to be applied to
 * every InvestmentBlock; while the -S option specifies a BlockSolverConfig
 * file for every InvestmentBlock. If the -B option is not provided when the
 * given netCDF file is a BlockFile, then a default configuration is
 * considered.
 *
 * Initial cuts can be provided by using the -l option. This option must be
 * followed by the path to the file containing the initial cuts. This file
 * must have the following format. The first line contains a header and its
 * content is ignored. Each of the following lines represent a cut and has the
 * following format:
 *
 *     t, a_0, a_1, ..., a_k, b
 *
 * where t is a stage (an integer between 0 and time horizon minus 1), a_0,
 * ..., a_k are the coefficients of the cut, and b is the constant term of the
 * cut.
 *
 * As a preprocessing, given redundant cuts can be removed by using the -e
 * option. Notice that all cuts will be subject to being removed, whether they
 * are provided in a netCDF file or by the -l option.
 *
 * There are a few ways to specify the initial state for the first stage
 * subproblem. This can be done by setting the initial state variable of
 * SDDPBlock or by setting the initial state parameter of SDDPGreedySolver.
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Rafael Durbano Lobato
 */

#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <queue>

#include <BatteryUnitBlock.h>
#include <BendersBlock.h>
#include <BlockSolverConfig.h>
#include <HydroSystemUnitBlock.h>
#include <IntermittentUnitBlock.h>
#include <NetworkBlock.h>
#include <RBlockConfig.h>
#include <SDDPBlock.h>
#include <StochasticBlock.h>
#include <SDDPGreedySolver.h>
#include <SDDPSolver.h>
#include <SlackUnitBlock.h>
#include <ThermalUnitBlock.h>
#include <UCBlock.h>

#include "CutProcessing.h"
#include "InvestmentBlock.h"
#include "InvestmentFunction.h"
#include "SDDPBlockSolutionOutput.h"

#ifdef USE_MPI
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#endif

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/

std::string filename{};
std::string block_config_filename{};
std::string solver_config_filename{};
std::string config_filename_prefix{};
std::string cuts_filename{};
std::string initial_point_filename{};

// State to be loaded into the InvestmentBlock Solver
std::string solver_state_input_filename{};

// Prefix to the name of the file that will store the State of the
// InvestmentBlock Solver
std::string solver_state_output_filename{};

const std::string best_solution_filename = "Solution_OUT.csv";

long num_sub_blocks_per_stage = 1;

bool eliminate_redundant_cuts = false;
bool simulate_investment = false;
bool single_scenario = false;
bool output_solution = false;
const bool force_hard_components = false;
const bool continuous_relaxation = true;

// Since BundleSolver cannot currently handle general bounds on the variables
// of the form l <= x <= u, these constraints must be reformulated by
// replacing them by 0 <= x <= u - l.
const bool reformulate_variable_bounds = true;

// This variable indicates whether negative prices may occur
const bool negative_prices = false;

std::string exe{};         ///< Name of the executable file
std::string docopt_desc{}; ///< Tool description

std::vector< double > initial_point;

/*--------------------------------------------------------------------------*/

// Gets the name of the executable from its full path
std::string get_filename( const std::string & fullpath ) {
 std::size_t found = fullpath.find_last_of( "/\\" );
 return( fullpath.substr( found + 1 ) );
}

/*--------------------------------------------------------------------------*/

void print_help() {
 // http://docopt.org
 std::cout << docopt_desc << std::endl;
 std::cout << "Usage:\n"
           << "  " << exe << " [options] <file>\n"
           << "  " << exe << " -h | --help\n"
           << std::endl
           << "Options:\n"
           << "  -a, --save-state <prefix>       Save states of the InvestmentBlock solver.\n"
           << "  -B, --blockcfg <file>           Block configuration.\n"
           << "  -b, --load-state <file>         Load a state for the InvestmentBlock solver.\n"
           << "  -c, --configdir <path>          The prefix for all config filenames.\n"
           << "  -e, --eliminate-redundant-cuts  Eliminate given redundant cuts.\n"
           << "  -h, --help                      Print this help.\n"
           << "  -l, --load-cuts <file>          Load cuts from a file.\n"
           << "  -n, --num-blocks <number>       Number of sub-Blocks per stage.\n"
           << "  -o, --output-solution           Output the solutions.\n"
           << "  -p, --prefix <path>             The prefix for all Block filenames.\n"
           << "  -S, --solvercfg <file>          Solver configuration.\n"
           << "  -s, --simulate                  Simulate the given investment.\n"
           << "  -x, --initial-investment <file> Initial investment."
           << std::endl;
}

/*--------------------------------------------------------------------------*/

long get_long_option() {
 char * end = nullptr;
 errno = 0;
 long option = std::strtol( optarg , &end , 10 );
 if( ( ! optarg ) || ( ( option = std::strtol( optarg , &end , 10 ) ) ,
                       ( errno || ( end && *end ) ) ) ) {
  option = -1;
 }
 return( option );
}

/*--------------------------------------------------------------------------*/

void process_args( int argc , char ** argv ) {

 if( argc < 2 ) {
  std::cout << exe << ": no input file\n"
            << "Try " << exe << "' --help' for more information.\n";
  exit( 1 );
 }

 const char * const short_opts = "a:B:b:c:hel:n:op:rS:sx:";
 const option long_opts[] = {
  { "save-state" ,               required_argument , nullptr , 'a' } ,
  { "blockcfg" ,                 required_argument , nullptr , 'B' } ,
  { "load-state" ,               required_argument , nullptr , 'b' } ,
  { "configdir" ,                required_argument , nullptr , 'c' } ,
  { "help" ,                     no_argument ,       nullptr , 'h' } ,
  { "eliminate-redundant-cuts" , no_argument ,       nullptr , 'e' } ,
  { "load-cuts" ,                required_argument , nullptr , 'l' } ,
  { "num-blocks" ,               required_argument , nullptr , 'n' } ,
  { "output-solution" ,          no_argument ,       nullptr , 'o' } ,
  { "prefix" ,                   required_argument , nullptr , 'p' } ,
  { "relax" ,                    no_argument ,       nullptr , 'r' } ,
  { "solvercfg" ,                required_argument , nullptr , 'S' } ,
  { "simulate" ,                 no_argument ,       nullptr , 's' } ,
  { "initial-investment" ,       required_argument , nullptr , 'x' } ,
  { nullptr ,                    no_argument ,       nullptr , 0 }
 };

 // Options
 while( true ) {
  const auto opt = getopt_long( argc , argv , short_opts ,
                                long_opts , nullptr );

  if( opt == -1 ) {
   break;
  }

  switch( opt ) {
   case 'a':
    solver_state_output_filename = std::string( optarg );
    break;
   case 'B':
    block_config_filename = std::string( optarg );
    break;
   case 'b':
    solver_state_input_filename = std::string( optarg );
    break;
   case 'c':
    config_filename_prefix = std::string( optarg );
    Configuration::set_filename_prefix( std::string( optarg ) );
    break;
   case 'e':
    eliminate_redundant_cuts = true;
    break;
   case 'l':
    cuts_filename = std::string( optarg );
    break;
   case 'n': {
    num_sub_blocks_per_stage = get_long_option();
    if( num_sub_blocks_per_stage <= 0 ) {
     std::cout << "The number of sub-Blocks per stage must be a "
               << "positive integer." << std::endl;
     exit( 1 );
    }
    break;
   }
   case 'o':
    output_solution = true;
    break;
   case 'p':
    Block::set_filename_prefix( std::string( optarg ) );
    break;
   case 'r':
    std::cout << "The -r option no longer exists. In order relax the "
              << "integrality constraints,\nplease properly configure the "
              << "solver. For instance, some solvers have the\nparameter "
              << "'intRelaxIntVars', which can be set to 1 in the solver\n"
              << "configuration file associated with the Block whose "
              << "constraints must be\nrelaxed." << std::endl;
    exit( 1 );
   case 'S':
    solver_config_filename = std::string( optarg );
    break;
   case 's': {
    simulate_investment = true;
    break;
   }
   case 'x':
    initial_point_filename = std::string( optarg );
    break;
   case 'h': // -h or --help
    print_help();
    exit( 0 );
   case '?': // Unrecognized option
   default:
    std::cout << "Try " << exe << "' --help' for more information.\n";
    exit( 1 );
  }
 }

 // Last argument
 if( optind < argc ) {
  filename = std::string( argv[ optind ] );
 }
 else {
  std::cout << exe << ": no input file\n"
            << "Try " << exe << "' --help' for more information.\n";
  exit( 1 );
 }
}

/*--------------------------------------------------------------------------*/

Block * get_uc_block( const SDDPBlock * sddp_block , const Index stage ) {

 auto benders_block = static_cast< BendersBlock * >
  ( sddp_block->get_sub_Block( stage )->get_inner_block() );

 auto objective = static_cast< FRealObjective * >
  ( benders_block->get_objective() );

 auto benders_function = static_cast< BendersBFunction * >
  ( objective->get_function() );

 return( benders_function->get_inner_block() );
}

/*--------------------------------------------------------------------------*/

bool update_hydro_unit( Block * previous_block , Block * block ,
                        const Index stage ) {
 auto unit = dynamic_cast< HydroUnitBlock * >( block );
 auto previous_unit = dynamic_cast< HydroUnitBlock * >( previous_block );

 if( ( ! unit ) && ( ! previous_unit ) )
  return( false );

 if( ( ! unit ) || ( ! previous_unit ) )
  throw( std::logic_error
         ( "sddp_solver: UCBlocks at stages " + std::to_string( stage - 1 ) +
           " and " + std::to_string( stage ) +
           " do not have the same structure." ) );

 auto number_generators = previous_unit->get_number_generators();

 if( number_generators != unit->get_number_generators() )
  throw( std::logic_error
         ( "sddp_solver: HydroUnitBlock at stage " +
           std::to_string( stage - 1 ) + " has " +
           std::to_string( number_generators ) +
           ", but corresponding HydroUnitBlock at stage " +
           std::to_string( stage ) + " has " +
           std::to_string( unit->get_number_generators() ) ) );

 const auto time_horizon = previous_unit->get_time_horizon();

 std::vector< double > flow_rate( number_generators );

 for( Index g = 0 ; g < number_generators ; ++g )
  flow_rate[ g ] =
   previous_unit->get_flow_rate( g , time_horizon - 1 )->get_value();

 unit->set_initial_flow_rate( flow_rate.cbegin() );

 return( true );
}

/*--------------------------------------------------------------------------*/

bool update_battery_unit( Block * previous_block , Block * block ,
                          const Index stage ) {
 auto unit = dynamic_cast< BatteryUnitBlock * >( block );
 auto previous_unit = dynamic_cast< BatteryUnitBlock * >( previous_block );

 if( ( ! unit ) && ( ! previous_unit ) )
  return( false );

 if( ( ! unit ) || ( ! previous_unit ) )
  throw( std::logic_error
         ( "sddp_solver: UCBlocks at stages " + std::to_string( stage - 1 ) +
           " and " + std::to_string( stage ) +
           " do not have the same structure." ) );

 const auto time_horizon = previous_unit->get_time_horizon();

 std::vector< double > initial_power_data =
  { ( previous_unit->get_active_power( 0 ) + time_horizon - 1 )->get_value() };

 unit->set_initial_power( initial_power_data.cbegin() );

 std::vector< double > initial_storage_data =
  { previous_unit->get_storage_level()[ time_horizon - 1 ].get_value() };

 unit->set_initial_storage( initial_storage_data.cbegin() );

 return( true );
}

/*--------------------------------------------------------------------------*/

int compute_init_up_down_time( const SDDPBlock * sddp_block ,
                               ThermalUnitBlock * previous_unit ,
                               ThermalUnitBlock * unit , const Index stage ) {

 auto time_horizon = previous_unit->get_time_horizon();
 auto commitment = previous_unit->get_commitment( 0 ) + time_horizon - 1;

 auto shutdown = previous_unit->get_shut_down( time_horizon - 1 );
 if( shutdown && shutdown->get_value() >= 0.5 ) {
  return( 0 );
 }

 int init_up_down_time = 0;
 const bool on = commitment->get_value() >= 0.5;
 if( on )
  init_up_down_time = 1;
 else
  init_up_down_time = -1;

 AbstractPath path;

 for( Index outer_t = 0 ; outer_t < stage ; ++outer_t ) {

  for( Index t = 1 ; t < time_horizon ; ++t , --commitment ) {
   if( std::abs( commitment->get_value() -
                 ( commitment - 1 )->get_value() ) > 0.5 )
    return( init_up_down_time );
   if( on ) ++init_up_down_time;
   else --init_up_down_time;
  }

  if( outer_t == stage - 1 )
   break;

  if( path.empty() ) {
   auto uc_block = get_uc_block( sddp_block , stage );
   path.build( unit , uc_block );
  }

  auto previous_uc_block = get_uc_block( sddp_block , stage - outer_t - 2 );
  previous_unit = dynamic_cast< ThermalUnitBlock * >
   ( path.get_element< Block >( previous_uc_block ) );

  time_horizon = previous_unit->get_time_horizon();

  if( ! previous_unit )
   throw( std::logic_error
          ( "sddp_solver::update_thermal_block: ThermalUnitBlock not found "
            "at stage " + std::to_string( stage - outer_t - 2 ) ) );

  commitment = previous_unit->get_commitment( 0 ) + time_horizon - 1;

  if( on ) {
   if( commitment->get_value() >= 0.5 ) ++init_up_down_time;
   else break;
  }
  else {
   if( commitment->get_value() < 0.5 ) --init_up_down_time;
   else break;
  }
 }

 return( init_up_down_time );
}

/*--------------------------------------------------------------------------*/

bool update_thermal_unit( const SDDPBlock * sddp_block ,
                          Block * previous_block , Block * block ,
                          const Index stage ) {

 auto previous_unit = dynamic_cast< ThermalUnitBlock * >( previous_block );
 auto unit = dynamic_cast< ThermalUnitBlock * >( block );

 if( ( ! unit ) && ( ! previous_unit ) )
  return( false );

 if( ( ! unit ) || ( ! previous_unit ) )
  throw( std::logic_error
         ( "sddp_solver: UCBlocks at stages " + std::to_string( stage - 1 ) +
           " and " + std::to_string( stage ) +
           " do not have the same structure." ) );

 if( simulate_investment && single_scenario ) {
  // One case where we can safely update the initial up and downtime is when
  // we are simulating a single scenario considering the simulation-based
  // investment function.
  auto init_up_down_time = compute_init_up_down_time
   ( sddp_block , previous_unit , unit , stage );

  std::vector< int > init_up_down_time_data = { init_up_down_time };
  unit->set_init_updown_time( init_up_down_time_data.cbegin() );
 }

 const auto time_horizon = previous_unit->get_time_horizon();

 std::vector< double > active_power_data =
  { ( previous_unit->get_active_power( 0 ) + time_horizon - 1 )->get_value() };
 unit->set_initial_power( active_power_data.cbegin() );

 return( true );
}

/*--------------------------------------------------------------------------*/

void callback( SDDPBlock * sddp_block , Block::Index stage ) {

 if( stage == 0 )
  return;

 auto previous_uc_block = get_uc_block( sddp_block , stage - 1 );
 auto uc_block = get_uc_block( sddp_block , stage );

 std::queue< Block * > blocks;
 blocks.push( uc_block );

 std::queue< Block * > previous_blocks;
 previous_blocks.push( previous_uc_block );

 while( ! blocks.empty() ) {
  auto block = blocks.front();
  blocks.pop();

  auto previous_block = previous_blocks.front();
  previous_blocks.pop();

  auto n = block->get_number_nested_Blocks();

  if( n != previous_block->get_number_nested_Blocks() ) {
   throw( std::logic_error
          ("sddp_solver: UCBlocks at stages " + std::to_string( stage - 1 ) +
           " and " + std::to_string( stage ) +
           " do not have the same structure." ) );
  }

  for( decltype( n ) i = 0 ; i < n ; ++i ) {
   blocks.push( block->get_nested_Block( i ) );
   previous_blocks.push( previous_block->get_nested_Block( i ) );
  }

  update_hydro_unit( previous_block , block , stage )
   || update_thermal_unit( sddp_block , previous_block , block , stage )
   || update_battery_unit( previous_block , block , stage );
 }
}

/*--------------------------------------------------------------------------*/

std::vector< double > get_default_initial_point( InvestmentBlock * block ) {
 block->generate_abstract_constraints();
 const auto & box_constraints = block->get_constraints();
 std::vector< double > initial_point( box_constraints.size() );
 for( Index i = 0 ; i < box_constraints.size() ; ++i ) {
  if( box_constraints[ i ].get_lhs() > -Inf< double >() )
   initial_point[ i ] = box_constraints[ i ].get_lhs();
  else if( box_constraints[ i ].get_rhs() < Inf< double >() )
   initial_point[ i ] = box_constraints[ i ].get_rhs();
  else
   initial_point[ i ] = 0;
 }

 return( initial_point );
}

/*--------------------------------------------------------------------------*/

std::vector< double > load_initial_point() {
 if( initial_point_filename.empty() )
  return{};

 std::ifstream file( initial_point_filename );

 // Make sure the file is open
 if( ! file.is_open() )
  throw( std::runtime_error( "It was not possible to open the file \"" +
                             initial_point_filename + "\"." ) );

 std::vector< double > initial_point;

 double component;
 while( file >> component )
  initial_point.push_back( component );

 return( initial_point );
}

/*--------------------------------------------------------------------------*/

void set_initial_point( InvestmentBlock * investment_block ) {

 // Generate the abstract variables so that we can set their values.

 investment_block->generate_abstract_variables();

 // Possibly load a given initial point.

 initial_point = load_initial_point();

 if( ! initial_point.empty() ) {
  // An initial point has been provided.

  const auto num_variables = investment_block->get_number_variables();
  if( initial_point.size() != num_variables )
   throw( std::logic_error( "The initial point has size " +
                            std::to_string( initial_point.size() ) + ", but "
                            "there are " + std::to_string( num_variables ) +
                            " variables." ) );

  if( reformulate_variable_bounds ) {

   // If variable bounds have been reformulated, the initial point must be
   // adjusted.

   const auto & var_lower_bound = investment_block->get_variable_lower_bound();
   for( Index i = 0 ; i < initial_point.size() ; ++i ) {
    if( ( i < var_lower_bound.size() ) &&
        ( var_lower_bound[ i ] > -Inf< double >() ) )
     initial_point[ i ] -= var_lower_bound[ i ];
   }
  }
 }
 else {
  // Since no initial point has been provided, we use the default one.
  initial_point = get_default_initial_point( investment_block );
 }

 // Finally, set the initial point.

 investment_block->set_variable_values( initial_point );
}

/*--------------------------------------------------------------------------*/

int get_objective_sense( const SDDPBlock * sddp_block ) {
 if( ! sddp_block )
  return( Objective::eUndef );
 return( sddp_block->get_objective_sense() );
}

/*--------------------------------------------------------------------------*/

void invest( InvestmentBlock * investment_block ) {

 auto investment_function = static_cast< InvestmentFunction * >
  ( investment_block->get_function() );

 // Possibly output the solution
 investment_function->
  set_par( InvestmentFunction::intOutputSolution , output_solution );

 for( Index i = 0 ; i < investment_function->get_number_nested_Blocks() ;
      ++i ) {

  auto sddp_block =
   dynamic_cast< SDDPBlock * >( investment_function->get_nested_Block( i ) );

  if( ! sddp_block ) {
   std::cout << "The sub-Block of the InvestmentBlock is not an SDDPBlock."
             << std::endl;
   exit( 1 );
  }

  auto sddp_solver = dynamic_cast< SDDPGreedySolver * >
   ( sddp_block->get_registered_solvers().front() );

  if( ! sddp_solver )
   throw( std::logic_error( "The Solver for the SDDPBlock must be an "
                            "SDDPGreedySolver." ) );

  sddp_solver->set_callback( [ sddp_block ]( Index stage ) {
   callback( sddp_block , stage );
  } );

  // Check whether there is a single scenario

  single_scenario = ( sddp_block->get_scenario_set().size() == 1 );

  // Load possibly given cuts

  if( ! cuts_filename.empty() ) {
   sddp_solver->set_par( SDDPGreedySolver::strLoadCuts , cuts_filename );
   sddp_solver->set_par( SDDPGreedySolver::intLoadCutsOnce , 1 );
  }

  // Eliminate redundant cuts if it is desired

  if( eliminate_redundant_cuts )
   CutProcessing().remove_redundant_cuts( sddp_block );
 }

 // Solve

 if( simulate_investment ) {

  if( ! initial_point.empty() ) {
   std::cout << "Simulating the investment (";

   const auto & var_lower_bound = investment_block->get_variable_lower_bound();

   for( Index i = 0 ; i < initial_point.size() ; ++i ) {

    auto x_i = initial_point[ i ];
    if( reformulate_variable_bounds && ( i < var_lower_bound.size() ) &&
        ( var_lower_bound[ i ] > -Inf< double >() ) )
     // Since variable bounds have been reformulated, adjust x_i so that the
     // user sees the expected initial point.
     x_i += var_lower_bound[ i ];

    if( i > 0 )
     std::cout << ", ";
    std::cout << x_i;
   }
   std::cout << ")." << std::endl;
  }

  // Disable the computation of linearization
  investment_function->
   set_par( InvestmentFunction::intComputeLinearization , 0 );

  // Simulate
  auto objective =
   static_cast< FRealObjective * >( investment_block->get_objective() );
  objective->compute();

  const auto value = objective->value();
  std::cout << "Value: " << std::setprecision( 20 ) << value << std::endl;
 }
 else {

  // Optimize

  auto investment_solver = investment_block->get_registered_solvers().front();

  investment_solver->set_log( &std::cout );

  // Output the variable and function values at each iteration
  investment_function->set_par( InvestmentFunction::strOutputFilename ,
                                "investment_candidates.txt" );

  if( ! solver_state_input_filename.empty() ) {
   // Load the given State

   netCDF::NcFile file;
   try {
    file.open( solver_state_input_filename , netCDF::NcFile::read );
    auto state = State::new_State( file );
    investment_solver->put_State( *state );
    delete( state );
   } catch( netCDF::exceptions::NcException & e ) {
    std::cout << "Warning: It was not possible to open the State file '"
              << solver_state_input_filename << "'. The State of the Solver "
              << "will not be loaded." << std::endl;
   } catch( const std::exception& e ) {
    std::cout << "Warning: An error occurred while loading the Solver State: '"
              << e.what() << "'." << std::endl;
   }
  }

  if( ! solver_state_output_filename.empty() ) {
   // Register an event to save the State of the Solver

   investment_solver->set_par( ThinComputeInterface::intEverykIt , 1 );
   investment_solver->set_event_handler
    ( ThinComputeInterface::eEverykIteration ,
      [ investment_solver ]() {
       static int i = 0;
       std::string filename =
        solver_state_output_filename + std::to_string( i++ ) + ".nc4";
       i %= 2;
       netCDF::NcFile file( filename , netCDF::NcFile::replace );
       investment_solver->serialize_State( file );
       return( ThinComputeInterface::eContinue );
      } );
  }

  // Register an event to keep track of the best solution and to handle the
  // output.

  std::vector< double > best_solution;

  const auto objective_sense =
   get_objective_sense( investment_function->get_sddp_block( 0 ) );

  const auto sign = ( objective_sense == Objective::eMin ) ? 1 : -1;
  const auto worst_value = sign * Inf< InvestmentFunction::FunctionValue >();
  auto best_solution_value = worst_value;

  investment_function->set_par( ThinComputeInterface::eBeforeTermination , 1 );
  investment_function->set_event_handler
   ( ThinComputeInterface::eBeforeTermination ,
     [ investment_block , investment_function , &best_solution_value ,
       &best_solution , objective_sense ]() {

      auto solution_improved = [ investment_function , &best_solution_value ,
                                 objective_sense ]() {
       const auto function_value = investment_function->get_value();
       if( ( ( objective_sense == Objective::eMin ) &&
             ( function_value < best_solution_value ) ) ||
           ( ( objective_sense == Objective::eMax ) &&
             ( function_value > best_solution_value ) ) ) {
        best_solution_value = function_value;
        return( true );
       }
       return( false );
      };

      if( solution_improved() ) {
       // Save the best solution found so far

       std::ofstream best_solution_file( best_solution_filename ,
                                         std::ios::out );

       const auto & variables = investment_block->get_variables();
       const auto & var_lower_bound =
        investment_block->get_variable_lower_bound();
       best_solution.resize( variables.size() );

       for( Index i = 0 ; i < variables.size() ; ++i ) {
        auto value = variables[ i ].get_value();
        if( reformulate_variable_bounds && ( i < var_lower_bound.size() ) &&
            ( var_lower_bound[ i ] > -Inf< double >() ) )
         value += var_lower_bound[ i ];
        best_solution[ i ] = value;
        best_solution_file << std::setprecision( 20 ) << value << std::endl;
       }

       // Possibly output information associated with the solution
       if( output_solution )
        SDDPBlockSolutionOutput().copy
         ( investment_function->get_sddp_block( 0 ) , ".best" , true );
      }

      return( ThinComputeInterface::eContinue );
     } );

  // Solve the investment problem

  investment_solver->compute();

#ifdef USE_MPI
  boost::mpi::communicator world;
  if( world.rank() == 0 ) {
#endif

   // Rename the output files if necessary

   if( output_solution )
    SDDPBlockSolutionOutput().rename
     ( investment_function->get_sddp_block( 0 ) , ".best" , true );

   // Output solution information

   if( best_solution_value == worst_value )
    std::cout << "No solution has been found." << std::endl;
   else if( best_solution_value == - worst_value )
    std::cout << "The problem is unbounded." << std::endl;
   else {
    // A solution has been found
    std::cout << "Solution value: " << std::setprecision( 20 )
              << best_solution_value << std::endl;
    const auto & var_lower_bound = investment_block->get_variable_lower_bound();
    const auto width = std::to_string( best_solution.size() ).size();
    for( Index i = 0 ; i < best_solution.size() ; ++i ) {
     auto value = best_solution[ i ];
     if( reformulate_variable_bounds && ( i < var_lower_bound.size() ) &&
         ( var_lower_bound[ i ] > -Inf< double >() ) )
      value += var_lower_bound[ i ];
     std::cout << std::setw( width ) << i << " " << value << std::endl;
    }
   }

#ifdef USE_MPI
  }
#endif
 }
}

/*--------------------------------------------------------------------------*/

void configure_Blocks( SDDPBlock * sddp_block ,
                       bool add_reserve_variables_to_objective ) {
 for( auto sub_block : sddp_block->get_nested_Blocks() ) {

  auto stochastic_block = static_cast< StochasticBlock * >( sub_block );
  auto benders_block = static_cast< BendersBlock * >
   ( stochastic_block-> get_nested_Blocks().front() );
  auto objective = static_cast< FRealObjective * >
   ( benders_block->get_objective() );
  auto benders_function = static_cast< BendersBFunction * >
   ( objective->get_function() );
  auto inner_block = benders_function->get_inner_block();

  std::queue< Block * > blocks;
  blocks.push( inner_block );

  while( ! blocks.empty() ) {
   auto block = blocks.front();
   blocks.pop();
   auto n = block->get_number_nested_Blocks();
   for( decltype( n ) i = 0 ; i < n ; ++i ) {
    blocks.push( block->get_nested_Block( i ) );
   }

   int cons_type = 1; // generate OneVarConstraints

   // Configure PolyhedralFunctionBlock
   if( auto polyhedral = dynamic_cast< PolyhedralFunctionBlock * >( block ) ) {
    auto config = new BlockConfig;
    config->f_static_variables_Configuration =
     new SimpleConfiguration< int >( 1 );
    polyhedral->set_BlockConfig( config );
   }

   else if( auto unit = dynamic_cast< SlackUnitBlock * >( block ) ) {
    auto config = new BlockConfig;
    config->f_static_constraints_Configuration =
     new SimpleConfiguration< int >( cons_type );
    unit->set_BlockConfig( config );
   }

   else if( auto unit = dynamic_cast< BatteryUnitBlock * >( block ) ) {
    auto config = new BlockConfig;
    config->f_static_variables_Configuration =
      new SimpleConfiguration< int >( negative_prices );
    config->f_static_constraints_Configuration =
     new SimpleConfiguration< int >( cons_type );
    unit->set_BlockConfig( config );
   }

   else if( auto unit = dynamic_cast< ThermalUnitBlock * >( block ) ) {
    auto config = new BlockConfig;
    config->f_static_constraints_Configuration =
     new SimpleConfiguration< int >( cons_type );

    if( add_reserve_variables_to_objective )
     config->f_objective_Configuration = new SimpleConfiguration< int >( 3 );

    unit->set_BlockConfig( config );
   }

  }
 }
}

/*--------------------------------------------------------------------------*/

void set_log( SDDPBlock * sddp_block , std::ostream * output_stream ) {
 for( auto sub_block : sddp_block->get_nested_Blocks() ) {
  auto stochastic_block = static_cast< StochasticBlock * >( sub_block );
  auto benders_block = static_cast< BendersBlock * >
   ( stochastic_block-> get_nested_Blocks().front() );
  auto objective = static_cast< FRealObjective * >
   ( benders_block->get_objective() );
  auto benders_function = static_cast< BendersBFunction * >
   ( objective->get_function() );
  auto inner_block = benders_function->get_inner_block();
  for( auto solver : inner_block->get_registered_solvers() )
   if( solver )
    solver->set_log( output_stream );
 }
}

/*--------------------------------------------------------------------------*/

void process_prob_file( const netCDF::NcFile & file ) {
 std::multimap< std::string , netCDF::NcGroup > problems = file.getGroups();
 // for each problem descriptor:
 for( auto & problem : problems ) {

  auto & problem_group = problem.second;

  // Deserialize the Block

  auto block_group = problem_group.getGroup( "Block" );
  auto block_type_att = block_group.getAtt( "type" );

  if( block_type_att.isNull() ) {
   std::cout << "The netCDF attribute 'type' was not found in the netCDF group "
             << block_group.getName() << "." << std::endl;
   exit( 1 );
  }

  std::string block_type;
  block_type_att.getValues( block_type );

  if( block_type != "InvestmentBlock" ) {
   std::cout << "The Block in the netCDF file " << block_type << " is "
             << block_type << ", but it must be an InvestmentBlock."
             << std::endl;
   exit( 1 );
  }

  auto investment_block = new InvestmentBlock;
  investment_block->set_number_sub_blocks( num_sub_blocks_per_stage );
  investment_block->deserialize( block_group );

  auto investment_function = static_cast< InvestmentFunction * >
   ( investment_block->get_function() );

  for( auto sddp_block_ : investment_function->get_nested_Blocks() ) {

   auto sddp_block = dynamic_cast< SDDPBlock * >( sddp_block_ );

   if( ! sddp_block ) {
    std::cout << "The sub-Block of the InvestmentBlock is not an SDDPBlock."
              << std::endl;
    exit( 1 );
   }

   // Set the output stream for the log of the inner Solvers

   set_log( sddp_block , &std::cout );

   // Eliminate redundant cuts if it is desired

   if( eliminate_redundant_cuts )
    CutProcessing().remove_redundant_cuts( sddp_block );
  }

  // Configure block
  auto block_config_group = problem_group.getGroup( "BlockConfig" );
  auto block_config = static_cast< BlockConfig * >
   ( BlockConfig::new_Configuration( block_config_group ) );
  if( ! block_config )
   throw( std::logic_error( "BlockConfig group was not properly provided." ) );
  block_config->apply( investment_block );
  block_config->clear();

  // Possibly set the initial point
  set_initial_point( investment_block );

  // Configure solver
  auto solver_config_group = problem_group.getGroup( "BlockSolver" );
  auto block_solver_config = static_cast< BlockSolverConfig * >
   ( BlockSolverConfig::new_Configuration( solver_config_group ) );
  if( ! block_solver_config )
   throw( std::logic_error( "BlockSolver group was not properly provided." ) );
  block_solver_config->apply( investment_block );
  block_solver_config->clear();

  std::cout << "Problem: " << problem.first << std::endl;

  // Solve

  invest( investment_block );

  // Destroy the Block and the Configurations

  block_config->apply( investment_block );
  delete( block_config );

  block_solver_config->apply( investment_block );
  delete( block_solver_config );

  delete( investment_block );
 }
}

/*--------------------------------------------------------------------------*/

BlockConfig * load_BlockConfig() {

 if( block_config_filename.empty() ) {
  std::cout << "Block configuration was not provided. "
   "Using default configuration." << std::endl;
  return( nullptr );
 }

 std::ifstream block_config_file;
 block_config_file.open( block_config_filename , std::ifstream::in );

 if( ! block_config_file.is_open() ) {
  std::cerr << "Block configuration " + block_config_filename +
   " was not found." << std::endl;
  exit( 1 );
 }

 std::cout << "Using Block configuration in " << block_config_filename
           << "." << std::endl;

 std::string config_name;
 block_config_file >> eatcomments >> config_name;
 auto config = Configuration::new_Configuration( config_name );
 auto block_config = dynamic_cast< BlockConfig * >( config );

 if( ! block_config ) {
  std::cerr << "Block configuration is not valid: "
            << config_name << std::endl;
  delete( config );
  exit( 1 );
 }

 try {
  block_config_file >> *block_config;
 }
 catch( const std::exception& e ) {
  std::cerr << "Block configuration is not valid: " << e.what() << std::endl;
  exit( 1 );
 }

 block_config_file.close();
 return( block_config );
}

/*--------------------------------------------------------------------------*/

BlockSolverConfig * load_BlockSolverConfig( const std::string & filename ) {

 if( filename.empty() ) {
  std::cout << "Solver configuration was not provided. "
   "Using default configuration." << std::endl;
  return( nullptr );
 }

 std::ifstream solver_config_file;
 solver_config_file.open( filename , std::ifstream::in );

 if( ! solver_config_file.is_open() ) {
  std::cerr << "Solver configuration " + filename +
   " was not found." << std::endl;
  exit( 1 );
 }

 std::cout << "Using Solver configuration in " << filename << "." << std::endl;

 std::string config_name;
 solver_config_file >> eatcomments >> config_name;
 auto config = Configuration::new_Configuration( config_name );
 auto solver_config = dynamic_cast< BlockSolverConfig * >( config );

 if( ! solver_config ) {
  std::cerr << "Solver configuration is not valid: "
            << config_name << std::endl;
  delete( config );
  exit( 1 );
 }

 try {
  solver_config_file >> *solver_config;
 }
 catch( ... ) {
  std::cout << "Solver configuration is not valid." << std::endl;
  exit( 1 );
 }

 solver_config_file.close();
 return( solver_config );
}

/*--------------------------------------------------------------------------*/

std::string get_str_par( ComputeConfig * compute_config ,
                         std::string par_name ) {
 for( const auto & pair : compute_config->str_pars ) {
  if( pair.first == par_name )
   return( pair.second );
 }
 return( "" );
}

/*--------------------------------------------------------------------------*/

int get_int_par( ComputeConfig * compute_config , std::string par_name ) {
 for( const auto & pair : compute_config->int_pars ) {
  if( pair.first == par_name )
   return( pair.second );
 }
 return( Inf< int >() );
}

/*--------------------------------------------------------------------------*/

bool using_lagrangian_dual_solver( BlockSolverConfig * sddp_solver_config ) {

 BlockSolverConfig * inner_solver_config = nullptr;
 ComputeConfig * compute_config = nullptr;

 for( Index i = 0 ; i < sddp_solver_config->num_ComputeConfig() ; ++i ) {

  if( sddp_solver_config->get_SolverName( i ) != "SDDPSolver" &&
      sddp_solver_config->get_SolverName( i ) != "ParallelSDDPSolver" &&
      sddp_solver_config->get_SolverName( i ) != "SDDPGreedySolver" )
   continue;

  compute_config = sddp_solver_config->get_SolverConfig( i );

  // Check if strInnerBSC is present

  auto strInnerBSC = get_str_par( compute_config , "strInnerBSC" );

  if( strInnerBSC.empty() )
   continue;

  // If it is, check if it is a config for a LagrangianDualSolver

  std::ifstream inner_solver_config_file
   ( config_filename_prefix + strInnerBSC , std::ifstream::in );

  if( ! inner_solver_config_file.is_open() )
   continue;

  std::string inner_config_name;
  inner_solver_config_file >> eatcomments >> inner_config_name;
  auto inner_config = Configuration::new_Configuration( inner_config_name );
  inner_solver_config = dynamic_cast< BlockSolverConfig * >( inner_config );

  if( ! inner_solver_config ) {
   inner_solver_config_file.close();
   delete( inner_config );
   continue;
  }

  try {
   inner_solver_config_file >> *inner_solver_config;
  }
  catch( ... ) {
   inner_solver_config_file.close();
   delete( inner_config );
   continue;
  }

  inner_solver_config_file.close();

  for( Index j = 0 ; j < inner_solver_config->num_ComputeConfig() ; ++j ) {
   if( inner_solver_config->get_SolverName( j ) == "LagrangianDualSolver" ) {
    delete( inner_config );
    return( true );
   }
  }
  delete( inner_config );
 }
 return( false );
}

/*--------------------------------------------------------------------------*/

void config_Lagrangian_dual( BlockSolverConfig * sddp_solver_config ,
                             SDDPBlock * sddp_block ,
                             InvestmentBlock * investment_block ) {

 if( sddp_block->get_number_nested_Blocks() == 0 )
  // The SDDPBlock has no sub-Block. There is nothing to be configured.
  return;

 BlockSolverConfig * inner_solver_config = nullptr;
 ComputeConfig * lagrangian_dual_compute_config = nullptr;
 ComputeConfig * compute_config = nullptr;

 // It indicates whether some Solver is a [Parallel]BundleSolver
 bool bundle_solver = false;
 bool do_easy_components = true;
 std::vector< int > vintNoEasy;

 // Index of the HydroSystemUnitBlock
 int hydro_system_index = -1;

 for( Index i = 0 ; i < sddp_solver_config->num_ComputeConfig() ; ++i ) {

  if( sddp_solver_config->get_SolverName( i ) != "SDDPSolver" &&
      sddp_solver_config->get_SolverName( i ) != "ParallelSDDPSolver" &&
      sddp_solver_config->get_SolverName( i ) != "SDDPGreedySolver" )
   continue;

  compute_config = sddp_solver_config->get_SolverConfig( i );

  // Check if strInnerBSC is present

  auto strInnerBSC = get_str_par( compute_config , "strInnerBSC" );

  if( strInnerBSC.empty() )
   return;

  // If it is, check if it is a config for a LagrangianDualSolver

  std::ifstream inner_solver_config_file
   ( config_filename_prefix + strInnerBSC , std::ifstream::in );

  if( ! inner_solver_config_file.is_open() )
   return;

  std::string inner_config_name;
  inner_solver_config_file >> eatcomments >> inner_config_name;
  auto inner_config = Configuration::new_Configuration( inner_config_name );
  inner_solver_config = dynamic_cast< BlockSolverConfig * >( inner_config );

  if( ! inner_solver_config ) {
   inner_solver_config_file.close();
   delete( inner_config );
   return;
  }

  try {
   inner_solver_config_file >> *inner_solver_config;
  }
  catch( ... ) {
   inner_solver_config_file.close();
   delete( inner_config );
   return;
  }

  inner_solver_config_file.close();

  for( Index j = 0 ; j < inner_solver_config->num_ComputeConfig() ; ++j ) {

   if( inner_solver_config->get_SolverName( j ) != "LagrangianDualSolver" )
    // It is not a ComputeConfig for a LagrangianDualSolver.
    // Check the next one.
    continue;

   lagrangian_dual_compute_config = inner_solver_config->get_SolverConfig( j );

   if( ! lagrangian_dual_compute_config )
    continue;

   // Find the inner Solver.
   auto sit = std::find_if( lagrangian_dual_compute_config->str_pars.begin() ,
                            lagrangian_dual_compute_config->str_pars.end() ,
                            []( auto & pair ) {
                             return( pair.first == "str_LDSlv_ISName" ); } );
   if( sit == lagrangian_dual_compute_config->str_pars.end() )
    // If it's not there, do nothing.
    continue;

   // Check if it is a [Parallel]BundleSolver.
   if( ( sit->second.find( "BundleSolver" ) == std::string::npos ) &&
       ( sit->second.find( "ParallelBundleSolver" ) == std::string::npos ) )
    continue;  // If it is not, do nothing.

   bundle_solver = true;

   // Check if the BundleSolver uses easy components.
   // Find if the ComputeConfig contains "intDoEasy".
   auto it = std::find_if( lagrangian_dual_compute_config->int_pars.begin() ,
                           lagrangian_dual_compute_config->int_pars.end() ,
                           []( auto & pair ) {
                            return( pair.first == "intDoEasy" ); } );
   if( it != lagrangian_dual_compute_config->int_pars.end() ) // if so
    do_easy_components = ( it->second & 1 ) > 0;  // read it
   else                               // otherwise
    do_easy_components = true;        // assume it is true (default)

   // We assume that there is at most one [Parallel]BundleSolver
   break;
  } // for each ComputeConfig for the inner Solver

  if( bundle_solver )
   break; // a BundleSolver has been found

 } // for each ComputeConfig for the Solver of SDDPBlock

 if( ! bundle_solver )
  // Since there is no BundleSolver, there is no need to configure any Block
  return;

 // The Configuration to be passed to get_var_solution() of the inner
 // Solver. We assume that only the HydroSystemBlock contains the necessary
 // part of the Solution (and that there is only one HydroSystemBlock) and
 // that the index of the HydroSystemBlock is the same at every stage.
 Configuration * get_var_solution_config = nullptr;

 // The Configuration to be passed to get_dual_solution() of the inner Solver.
 Configuration * get_dual_solution_config = nullptr;

 const std::string thermal_config_filename = "TUBSCfg.txt";
 const std::string hydro_config_filename = "HSUBSCfg.txt";
 const std::string other_unit_config_filename = "OUBSCfg.txt";
 const std::string default_config_filename = "LPBSCfg.txt";

 enum ConfigIndex { thermal = 0 , hydro , other_unit , default_config };

 // Vector with unique names of Configuration files ordered according to the
 // ConfigIndex enum.
 const std::vector< std::string > vstr_LDSl_Cfg = { thermal_config_filename ,
  hydro_config_filename , other_unit_config_filename ,
  default_config_filename };

 // We assume that all sub-Blocks of SDDPBlock have the same structure.

 const auto sub_block = sddp_block->get_nested_Block( 0 );

 auto stochastic_block = static_cast< StochasticBlock * >( sub_block );
 auto benders_block = static_cast< BendersBlock * >
  ( stochastic_block-> get_nested_Blocks().front() );
 auto objective = static_cast< FRealObjective * >
  ( benders_block->get_objective() );
 auto benders_function = static_cast< BendersBFunction * >
  ( objective->get_function() );
 auto inner_block = benders_function->get_inner_block();

 std::vector< int > vint_LDSl_WBSCfg;
 vint_LDSl_WBSCfg.reserve( inner_block->get_number_nested_Blocks() );

 /* The vector "required_primal_solution" will store the indices of Blocks
  * whose primal solutions are required (during the solution process). In
  * SDDP, only the primal solution of the HydroSystemUnitBlock is necessary
  * (as only the final volumes of the reservoirs are required during the
  * solution process). In simulation mode, the primal solutions that are
  * required are those of the Blocks that link two consecutive stages, which
  * are HydroSystemUnitBlock, BatteryUnitBlock, and ThermalUnitBlock.
  *
  * Notice that, in simulation mode, not all Blocks have their primal
  * solutions retrieved, which impacts the part of the solution that is output
  * (see UCBlockSolutionOutput). If the solutions of other Blocks are required
  * to be output when using LagrangianDualSolver+BundleSolver, then the
  * indices of these Blocks must be added to the vector
  * "required_primal_solution".
  *
  * This is currently not done due to a limitation of BundleSolver.
  * BundleSolver does not currently provide primal solutions for easy
  * components. Therefore, in order to have the primal solution of Blocks
  * other than HydroSystemUnitBlock, BatteryUnitBlock, and ThermalUnitBlock,
  * these Blocks must be treated as hard components (and they are currently
  * treated as easy components). Once BundleSolver is capable of providing
  * primal solutions of easy components, these Blocks can remain as easy
  * components and their indices can simply be added to the vector
  * "required_primal_solution". */

 std::vector< int > required_primal_solution;

 int inner_sub_block_index = 0;
 for( auto inner_sub_block : inner_block->get_nested_Blocks() ) {

  if( dynamic_cast< BatteryUnitBlock * >( inner_sub_block ) ) {

   required_primal_solution.push_back( inner_sub_block_index );

   // The primal solution of the BatteryUnitBlock is required as the storage
   // levels link two consecutive stages. Since BundleSolver currently does
   // not provide primal solutions for easy components, the BatteryUnitBlock
   // must be treated as a hard component. Once this feature is implemented by
   // BundleSolver, the BatteryUnitBlock can become an easy component.
   vint_LDSl_WBSCfg.push_back( ConfigIndex::other_unit );
   vintNoEasy.push_back( inner_sub_block_index );
  }
  if( dynamic_cast< ThermalUnitBlock * >( inner_sub_block ) ) {

   required_primal_solution.push_back( inner_sub_block_index );

   // ThermalUnitBlock is a non-easy component since there is a specialized
   // solver for it.
   vint_LDSl_WBSCfg.push_back( ConfigIndex::thermal );
   vintNoEasy.push_back( inner_sub_block_index );
  }
  else if( dynamic_cast< HydroSystemUnitBlock * >( inner_sub_block ) ) {
   required_primal_solution.push_back( inner_sub_block_index );
   hydro_system_index = inner_sub_block_index;

   // The HydroSystemUnitBlock could be treated as an easy component. However,
   // due to a current limitation of BundleSolver, the HydroSystemUnitBlock is
   // considered a hard component. This is because its primal solution (the
   // volume of the reservoirs) is required both in SDDP and in simulation
   // mode, but BundleSolver cannot currently provide primal solutions for
   // easy components. Once this feature is implemented by BundleSolver, the
   // HydroSystemUnitBlock can become an easy component.

   vint_LDSl_WBSCfg.push_back( ConfigIndex::hydro );
   vintNoEasy.push_back( inner_sub_block_index );
  }
  else if( dynamic_cast< IntermittentUnitBlock * >( inner_sub_block ) ) {
   required_primal_solution.push_back( inner_sub_block_index );
   vint_LDSl_WBSCfg.push_back( ConfigIndex::other_unit );
   vintNoEasy.push_back( inner_sub_block_index );
  }
  else if( dynamic_cast< NetworkBlock * >( inner_sub_block ) ) {
   /* TODO Dual solutions of the NetworkBlocks are necessary only if there are
    * transmission lines that are subject to investment. Since BundleSolver
    * currently does not provide solutions for easy components, the
    * NetworkBlock must be treated as a hard component. Once this feature is
    * implemented by BundleSolver, the NetworkBlock can become an easy
    * component. */
   vint_LDSl_WBSCfg.push_back( ConfigIndex::default_config );
   vintNoEasy.push_back( inner_sub_block_index );
  }
  else if( ! do_easy_components ) {
   vintNoEasy.push_back( inner_sub_block_index );
   if( dynamic_cast< UnitBlock * >( inner_sub_block ) )
    vint_LDSl_WBSCfg.push_back( ConfigIndex::other_unit );
   else
    vint_LDSl_WBSCfg.push_back( ConfigIndex::default_config );
  }
  else
   vint_LDSl_WBSCfg.push_back( ConfigIndex::default_config );

  ++inner_sub_block_index;
 }

 if( ! vintNoEasy.empty() ) {
  // Remove any vintNoEasy parameter that is possibly there
  lagrangian_dual_compute_config->vint_pars.erase
   ( std::remove_if( lagrangian_dual_compute_config->vint_pars.begin() ,
                     lagrangian_dual_compute_config->vint_pars.end() ,
                     []( const auto & pair ) {
                      return( pair.first == "vintNoEasy" ); } ) ,
     lagrangian_dual_compute_config->vint_pars.end() );

  // Add the vintNoEasy parameter that was constructed here
  lagrangian_dual_compute_config->vint_pars.push_back
   ( std::make_pair( "vintNoEasy" , std::move( vintNoEasy ) ) );
 }

 lagrangian_dual_compute_config->vint_pars.push_back
  ( std::make_pair( "vint_LDSl_WBSCfg" , std::move( vint_LDSl_WBSCfg ) ) );

 lagrangian_dual_compute_config->vstr_pars.push_back
  ( std::make_pair( "vstr_LDSl_Cfg" , std::move( vstr_LDSl_Cfg ) ) );

 // Configuration for the sub-Blocks may need to be cloned since the same
 // Configuration is used to configure multiple Blocks.
 lagrangian_dual_compute_config->int_pars.push_back
  ( std::make_pair( "int_LDSlv_CloneCfg" , 1 ) );

 compute_config->str_pars.erase
  ( std::remove_if( compute_config->str_pars.begin() ,
                    compute_config->str_pars.end() ,
                    []( const auto & pair ) {
                     return( pair.first == "strInnerBSC" ); } ) ,
    compute_config->str_pars.end() );

 /* The extra Configuration of the SDDPSolver and the SDDPGreedySolver is a
  * vector with pointers to the following elements (in that order):
  *
  * - a BlockConfig (which is currently nullptr) for the inner Block;
  *
  * - a BlockSolverConfig for the inner Block;
  *
  * - the Configuration to be passed to get_var_solution() when retrieving
  *   the Solutions to the inner Blocks of the BendersBFunctions.
  *
  * The extra Configuration of the SDDPGreedySolver has an additional (fourth)
  * element, which is
  *
  * - the Configuration to be passed to get_dual_solution() when retrieving
  *   the dual Solutions to the inner Blocks of the BendersBFunctions. */

 Configuration * extra_config = nullptr;

 /* Here we create a Configuration for
  * LagrangianDualSolver::get_var_solution() that requires the primal
  * solutions only of certain Blocks. In SDDP, only the primal solution of the
  * HydroSystemUnitBlock is necessary (as only the final volumes of the
  * reservoirs are required during the solution process). In simulation mode,
  * the solutions that are required are those of the Blocks that link two
  * consecutive stages, which are HydroSystemUnitBlock, BatteryUnitBlock, and
  * ThermalUnitBlock. */

 get_var_solution_config = new SimpleConfiguration< std::vector< int > >
  ( required_primal_solution );

 /* In investment mode, the only part of the dual solution that is required
  * is that associated with the UnitBlocks that are subject to
  * investment. Moreover, if transmission lines are also subject to
  * investment, then the dual solutions of all NetworkBlocks are also
  * necessary. */

 const auto ucblock = dynamic_cast< const UCBlock * >( inner_block );
 if( ! ucblock ) {
  std::cout << "In investment mode, the sub-problem must be a UCBlock."
            << std::endl;
  exit( 1 );
 }
 const auto time_horizon = ucblock->get_time_horizon();

 auto investment_function = static_cast< InvestmentFunction * >
  ( investment_block->get_function() );

 // Indices of the assets that are subject to investment.
 const auto & asset_indices = investment_function->get_asset_indices();

 // Types of assets that are subject to investment.
 const auto & asset_type = investment_function->get_asset_type();

 // Number of UnitBlocks and lines that are subject to investment.
 Index num_blocks = 0;
 Index num_lines = 0;
 for( const auto & type : asset_type ) {
  if( type == InvestmentFunction::eUnitBlock )
   ++num_blocks;
  else if( type == InvestmentFunction::eLine )
   ++num_lines;
 }

 // List containing the indices of the sub-Blocks of the UCBlock that are
 // subject to investment.
 std::vector< std::pair< int , int > > required_dual_solution;
 required_dual_solution.reserve( num_blocks + time_horizon );

 for( Index i = 0 ; i < asset_type.size() ; ++i )
  if( asset_type[ i ] == InvestmentFunction::eUnitBlock )
   required_dual_solution.push_back( { asset_indices[ i ] , -1 } );

 if( num_lines > 0 ) {
  // Since there are lines which are subject to investment, we must require
  // the dual solutions of all NetworkBlocks.

  const auto num_ucblock_sub_blocks = ucblock->get_number_nested_Blocks();
  Index num_network_blocks = 0;

  for( Index i = 0 ; i < num_ucblock_sub_blocks ; ++i ) {
   if( dynamic_cast< NetworkBlock * >( ucblock->get_nested_Block( i ) ) ) {
    required_dual_solution.push_back( { i , -1 } );
    ++num_network_blocks;
   }
  }

  // Check whether the number of NetworkBlocks is equal to the time horizon.

  if( num_network_blocks != time_horizon ) {
   std::cout << "The number of expected NetworkBlocks in the UCBlock is "
             << time_horizon << ", but " << num_network_blocks
             << " were found." << std::endl;
   exit( 1 );
  }
 } // end( if( num_lines > 0 ) )

 // To require the dual solution of the linking constraints, we add the pair
 // (-1, -1).
 required_dual_solution.push_back( { -1 , -1 } );

 // Finally create the SimpleConfiguration for the get_dual_solution() method.

 get_dual_solution_config = new SimpleConfiguration
  < std::vector< std::pair< int , int > > >( required_dual_solution );

 // Create the extra Configuration for SDDPGreedySolver.

 extra_config = new SimpleConfiguration< std::vector< Configuration * > >
  ( { nullptr , inner_solver_config , get_var_solution_config ,
     get_dual_solution_config } );

 // Set the extra Configuration

 compute_config->f_extra_Configuration = extra_config;

 // OSIMPSolver is currently not able to deal with some changes in a Block
 // (for instance, when some bound structure changes). In order to try to
 // avoid this case, we set a scenario, so that when OSIMPSolver is attached
 // to a Block, the data in that Block is a relevant one and, hopefully, will
 // not later be responsible for any other change in the bound structure. If
 // OSIMPSolver still complains, then other actions may be required (for
 // instance, replacing zeros by very small numbers in the scenarios).

 for( Index t = 0 ; t < sddp_block->get_time_horizon() ; ++t )
  for( Index i = 0 ; i < sddp_block->get_num_sub_blocks_per_stage() ; ++i )
   sddp_block->set_scenario( 0 , t , i );
}

/*--------------------------------------------------------------------------*/

void process_block_file( const netCDF::NcFile & file ) {
 std::multimap< std::string , netCDF::NcGroup > blocks = file.getGroups();

 // BlockConfig
 auto given_block_config = load_BlockConfig();

 BlockConfig * block_config = nullptr;
 if( given_block_config ) {
  block_config = given_block_config->clone();
  block_config->clear();
 }

 // BlockSolverConfig
 auto solver_config = load_BlockSolverConfig( solver_config_filename );
 if( ! solver_config ) {
  std::cout << "The Solver configuration is not valid." << std::endl;
  exit( 1 );
 }

 auto cleared_solver_config = solver_config->clone();
 cleared_solver_config->clear();

 // For each Block descriptor
 for( auto block_description : blocks ) {

  // Deserialize the Block

  auto block_type_att = block_description.second.getAtt( "type" );

  if( block_type_att.isNull() ) {
   std::cout << "The netCDF attribute 'type' was not found in the netCDF "
             << "group " << block_description.second.getName() << "." << std::endl;
   exit( 1 );
  }

  std::string block_type;
  block_type_att.getValues( block_type );

  if( block_type != "InvestmentBlock" ) {
   std::cout << "The Block in the netCDF file " << block_type << " is "
             << block_type << ", but it must be an InvestmentBlock."
             << std::endl;
   exit( 1 );
  }

  auto investment_block = new InvestmentBlock;
  investment_block->set_number_sub_blocks( num_sub_blocks_per_stage );
  investment_block->deserialize( block_description.second );

  auto investment_function = static_cast< InvestmentFunction * >
   ( investment_block->get_function() );

  for( auto sddp_block_ : investment_function->get_nested_Blocks() ) {

   auto sddp_block = dynamic_cast< SDDPBlock * >( sddp_block_ );

   if( ! sddp_block ) {
    std::cout << "The sub-Block of the InvestmentBlock is not an SDDPBlock."
              << std::endl;
    exit( 1 );
   }

   // Set the output stream for the log of the inner Solvers

   set_log( sddp_block , &std::cout );

   // Eliminate redundant cuts if it is desired

   if( eliminate_redundant_cuts )
    CutProcessing().remove_redundant_cuts( sddp_block );
  }

  // Configure the SDDPBlock

  bool is_using_lagrangian_dual_solver = false;

  if( given_block_config )
   given_block_config->apply( investment_block );
  else {
   is_using_lagrangian_dual_solver =
    using_lagrangian_dual_solver( solver_config );

   for( auto sddp_block_ : investment_function->get_nested_Blocks() ) {

    auto sddp_block = dynamic_cast< SDDPBlock * >( sddp_block_ );

    configure_Blocks( sddp_block , is_using_lagrangian_dual_solver );
   }

   if( reformulate_variable_bounds ) {
    // Since BundleSolver cannot currently handle general bounds on the
    // variables of the form l <= x <= u, we create a BlockConfig to instruct
    // the InvestmentBlock to reformulate the bound constraints by replacing
    // l <= x <= u by 0 <= x <= u - l.
    auto config = new BlockConfig;
    config->f_static_constraints_Configuration =
     new SimpleConfiguration< int >( 1 );

    investment_block->set_BlockConfig( config );
   }
  }

  // Configure the Solver

  if( is_using_lagrangian_dual_solver ) {
   auto investment_function = static_cast< InvestmentFunction * >
    ( investment_block->get_function() );
   for( auto sddp_block_ : investment_function->get_nested_Blocks() ) {
    auto sddp_block = dynamic_cast< SDDPBlock * >( sddp_block_ );
    config_Lagrangian_dual( solver_config , sddp_block , investment_block );
   }
  }

  // TODO This config file must be indicated in some appropriate way.

  const auto filename = config_filename_prefix + "sddp_greedy_investment.txt";

  auto sddp_solver_config = load_BlockSolverConfig( filename );

  ComputeConfig investment_function_config;

  investment_function_config.f_extra_Configuration =
   new SimpleConfiguration< std::map< std::string , Configuration * > >
   ( { { "BlockSolverConfig" , sddp_solver_config  } } );

  investment_function->set_ComputeConfig( &investment_function_config );

  // Possibly set the initial point

  set_initial_point( investment_block );

  // Finally, apply the Solver configuration

  solver_config->apply( investment_block );

  // Solve

  invest( investment_block );

  // Destroy the InvestmentBlock and the Configurations

  if( block_config )
   block_config->apply( investment_block );
  if( ! given_block_config ) {
   delete( block_config );
   block_config = nullptr;
  }

  cleared_solver_config->apply( investment_block );

  delete( investment_block );
 }

 delete( block_config );
 delete( given_block_config );
 delete( solver_config );
 delete( cleared_solver_config );
}

/*--------------------------------------------------------------------------*/

/// returns the final state (solution) of the system at the given \p stage
std::vector< double > get_final_state( SDDPBlock * block , Index stage ) {

 Index state_size = 0;
 for( Index i = 0 ; i < block->get_num_polyhedral_function_per_sub_block() ;
      ++i ) {
  state_size +=
   block->get_polyhedral_function( stage , i )->get_num_active_var();
 }

 std::vector< double > state;
 state.reserve( state_size );

 for( Index i = 0 ; i < block->get_num_polyhedral_function_per_sub_block() ;
      ++i ) {
  const auto polyhedral_function = block->get_polyhedral_function( stage , i );
  for( const auto & variable : * polyhedral_function ) {
   state.push_back
    ( static_cast< const ColVariable & >( variable ).get_value() );
  }
 }
 return( state );
}

/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv ) {

#ifdef USE_MPI
 boost::mpi::environment env( argc , argv );
#endif

 docopt_desc = "SMS++ investment solver.\n";
 exe = get_filename( argv[ 0 ] );
 process_args( argc , argv );

 if( solver_config_filename.empty() ) {
  // For the moment, the Solver configuration must be provided.
  std::cout << "A Solver configuration must be provided." << std::endl;
  exit( 0 );
 }

 netCDF::NcFile file;
 try {
  file.open( filename , netCDF::NcFile::read );
 } catch( netCDF::exceptions::NcException & e ) {
  std::cerr << "Cannot open nc4 file " << filename << std::endl;
  exit( 1 );
 }

 netCDF::NcGroupAtt gtype = file.getAtt( "SMS++_file_type" );
 if( gtype.isNull() ) {
  std::cerr << filename << " is not an SMS++ nc4 file." << std::endl;
  exit( 1 );
 }

 int type;
 gtype.getValues( &type );

 switch( type ) {
  case eProbFile: {
   std::cout << filename << " is a problem file, "
    "ignoring Block/Solver configurations..." << std::endl;
   process_prob_file( file );
   break;
  }

  case eBlockFile: {
   std::cout << filename << " is a block file." << std::endl;
   process_block_file( file );
   break;
  }

  default:
   std::cerr << filename << " is not a valid SMS++ file." << std::endl;
   exit( 1 );
 }

 file.close();
 return( 0 );
}
