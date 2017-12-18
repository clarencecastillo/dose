'''
File containing support functions for running a simulation. The functions 
in this file is not for public use; all functions in this file are private 
functions.

Date created: 10th October 2013
'''
import random, inspect, os
import os.path
from datetime import datetime
from time import time
from copy import deepcopy
from shutil import copyfile

# In Python 3, cPickle is no longer needed: Py3 looks for
# an optimized version, and if it founds none, will load the
# pure python implementation of pickle. 
try:
    import cPickle as pickle
except ImportError:
    import pickle

import dose_world
import genetic
import ragaraja, register_machine

from database_calls import connect_database, db_log_simulation_parameters
from database_calls import db_report

from . import events


def simulation_core(sim_functions, sim_parameters, eb, populations, world):
    '''
    Sequential ecological cell DOSE simulator.
    
    Performs the following operations:
        - Creating simulation file directory for results text file, population 
        freeze, and world burial storage
        - Generate a simulation start time to identify the current simulation
        - Define active Ragaraja instructions
        - Connecting to logging database (if needed)
        - Writing simulation parameters into results text file
        - Initialize World and Population
        - Deploy population(s) onto the world
        - Run the simulation and recording the results
        - Writing final results into results text file
        - Close logging database (if used)
        - Copy simulation script into simulation file directory

    Events available for subscription and data:
        - SIMULATION_START:
            - time_start: datetime = simulation start time
            - directory: String = simulation file path directory
            - max_generations: int = maximum generations
            - generation_count: int = starting generation
        - GENERATION_WORLD_UPDATE:
            - world: dose_world.World = world object
            - generation_count: int = current generation
        - GENERATION_POPULATIONS_UPDATE:
            - populations: dict = dictionary of population objects
            - generation_count: int = current generation
        - SIMULATION_END:
            - time_end: datetime = simulation end time
    
    @param sim_functions: implemented simulation functions (see 
    dose.dose_functions)
    @param sim_parameters: simulation parameters dictionary (see Examples)
    @param eb: An instance of dose.event_broker
    @param populations: dictionary of population objects
    @param world: dose_world.World object
    '''

    eb.log('Started simulation preparation')

    time_start = datetime.utcnow()
    time_start_str = '-'.join([str(time_start).split(' ')[0], str(time())])

    eb.log('Creating simulation file directories')
    directory = '_'.join([sim_parameters["simulation_name"], time_start_str])
    directory = os.sep.join([os.getcwd(), 'Simulations', directory]) + os.sep
    if not os.path.exists(directory):
        os.makedirs(directory)

    eb.log('Adding simulation directory to simulation parameters')
    sim_parameters["directory"] = directory

    eb.log('Adding starting time to simulation parameters')
    sim_parameters["starting_time"] = time_start_str

    sim_functions = sim_functions()

    eb.log('Activating ragaraja version: ' + str(sim_parameters["ragaraja_version"]))
    if sim_parameters["ragaraja_version"] == 0:
        ragaraja.activate_version(sim_parameters["ragaraja_version"], sim_parameters["ragaraja_instructions"])
    else:
        ragaraja.activate_version(sim_parameters["ragaraja_version"])

    con = None
    cur = None
    if "database_file" in sim_parameters and "database_logging_frequency" in sim_parameters:
        eb.log('Connecting to database file: ' + sim_parameters["database_file"])
        (con, cur) = connect_database(None, sim_parameters)
        eb.log('Logging simulation parameters to database file')
        (con, cur) = db_log_simulation_parameters(con, cur, sim_parameters)

    max_generations = None
    generation_count = None
    for pop_name in populations:
        eb.log('Preparing population: ' + pop_name + ' for simulation')
        if 'sim_folder' in sim_parameters or 'database_source' in sim_parameters:
            eb.log('Calculating final generation count')
            max_generations = sim_parameters["rev_start"][sim_parameters["population_names"].index(pop_name)] + \
                sim_parameters["extend_gen"]
            eb.log('Updating generation count from previous simulation')
            generation_count = sim_parameters["rev_start"][0]
        else:
            eb.log('Deploying population in World entity')
            if sim_parameters["deployment_code"] == 0:
                eb.log('Executing user defined deployment scheme')
                deploy_0(sim_parameters, populations, pop_name, world)
            elif sim_parameters["deployment_code"] == 1:
                eb.log('Executing deployment code 1: Single eco-cell deployment')
                deploy_1(sim_parameters, populations, pop_name, world)
            elif sim_parameters["deployment_code"] == 2:
                eb.log('Executing deployment code 2: Random eco-cell deployment')
                deploy_2(sim_parameters, populations, pop_name, world)
            elif sim_parameters["deployment_code"] == 3:
                eb.log('Executing deployment code 3: Even eco-cell deployment')
                deploy_3(sim_parameters, populations, pop_name, world)
            elif sim_parameters["deployment_code"] == 4:
                eb.log('Executing deployment code 4: Centralized eco-cell deployment')
                deploy_4(sim_parameters, populations, pop_name, world)
            max_generations = sim_parameters["maximum_generations"]
            generation_count = 0

        eb.log('Writing simulation parameters into txt file report')
        write_parameters(sim_parameters, pop_name)

        eb.log('Updating generation count')
        populations[pop_name].generation = generation_count

    eb.log('Simulation preparation complete')
    eb.publish(events.SIMULATION_START, {
        'time_start': time_start,
        'directory': directory,
        'max_generations': max_generations,
        'generation_count': generation_count
    })

    while generation_count < max_generations:
        generation_count += 1

        sim_functions.ecoregulate(world)
        eco_cell_iterator(world, sim_parameters, sim_functions.update_ecology)
        eco_cell_iterator(world, sim_parameters, sim_functions.update_local)
        eco_cell_iterator(world, sim_parameters, sim_functions.report)
        eb.publish(events.GENERATION_WORLD_UPDATE, {
            'world': world,
            'generation': generation_count
        })

        if generation_count % int(sim_parameters["eco_buried_frequency"]) == 0:
            bury_world(sim_parameters, world, generation_count)

        for pop_name in populations:
            if sim_parameters["interpret_chromosome"]:
                interpret_chromosome(sim_parameters, populations, pop_name, world)

            report_generation(sim_parameters, populations, pop_name, sim_functions, generation_count)
            sim_functions.organism_movement(populations, pop_name, world)
            sim_functions.organism_location(populations, pop_name, world)

        eb.publish(events.GENERATION_POPULATIONS_UPDATE, {
            'populations': populations,
            'generation': generation_count
        })

        if "database_file" in sim_parameters and "database_logging_frequency" in sim_parameters and \
                generation_count % int(sim_parameters["database_logging_frequency"]) == 0:
                (con, cur) = db_report(con, cur, sim_functions, sim_parameters["starting_time"], populations, world,
                                       generation_count)

        eb.log('Generation ' + str(generation_count) + ' complete')

    eb.log('Closing simulation results')
    for pop_name in populations:
        close_results(sim_parameters, pop_name)

    if "database_file" in sim_parameters and "database_logging_frequency" in sim_parameters:
        eb.log('Committing logged data into database file and terminating database connection')
        con.commit()
        con.close()

    eb.log('Copying simulation file script to simulation results directory')

    # DV on my system, the inspect.stack()[2][1] value returns the full
    # path to the file ('/home/douwe/scr/.../scriptname.py').
    # The end result is an error for the original code, below, as the
    # copyfile command is pointed to a wrong location: the directory
    # is added twice. This might be Py3, or Windows functionality.
    # With using the os.path module, this should give an OS-independent
    # way of extracting the basename and linking to the directory.
    sim_script_basename = os.path.basename(inspect.stack()[2][1])
    copyfile(inspect.stack()[2][1],
            # DV_original
            # sim_parameters['directory'] + inspect.stack()[2][1])
             os.path.join(sim_parameters['directory'], sim_script_basename))
    eb.log('Simulation ended')
    eb.publish(events.SIMULATION_END, {
        'time_end': datetime.utcnow()
    })


def coordinates(location):
    '''
    Helper function to transpose ecological cell into a tuple.
    
    @param location: location of ecological cell as 3-element iterable 
    data type
    @return: location of ecological cell as (x,y,z)
    '''

    x = location[0]
    y = location[1]
    z = location[2]
    return (x,y,z)


def adjacent_cells(sim_parameters, location):
    '''
    Function to get a list of adjacent ecological cells from a given 
    location.
    
    @param sim_parameters: simulation parameters dictionary (see Examples)
    @param location: location of ecological cell as (x,y,z)
    @return: list of locations of adjacent cells.
    '''
    trashbin = []
    temp_cells = []
    world_size = [sim_parameters["world_x"],
                  sim_parameters["world_y"],
                  sim_parameters["world_z"]]
    for i in range(3):
        new_location = [spot for spot in location]
        new_location[0] += 1
        temp_cells.append(new_location)
        new_location = [spot for spot in location]
        new_location[0] -= 1
        temp_cells.append(new_location)    
    for i in range(2):
        new_location = [spot for spot in location]
        new_location[1] -= 1
        temp_cells.append(new_location) 
    for i in range(0,4,3):
        temp_cells[i][1] += 1
        temp_cells[i+1][1] -= 1
    temp_cells[-1][1] += 2
    for i in range(8):
        for x in range(2):
            if temp_cells[i][x] >= world_size[x] or temp_cells[i][x] < 0:
                if temp_cells[i] not in trashbin:
                    trashbin.append(temp_cells[i])
    for location in trashbin:
        temp_cells.remove(location)
    return [tuple(location) for location in temp_cells]


def spawn_populations(sim_parameters):
    '''
    Initializing starting population(s) for a simulation. Each organism 
    in each population at this stage will be genetic clones of each 
    other.
    
    @param sim_parameters: simulation parameters dictionary (see Examples)
    @return: dictionary of population objects with population name as key
    '''
    temp_Populations = {}
    for pop_name in sim_parameters["population_names"]:
        temp_Populations[pop_name] = genetic.population_constructor(sim_parameters)
        for individual in temp_Populations[pop_name].agents:
            individual.generate_name()
            individual.status['deme'] = pop_name
    return temp_Populations


def eco_cell_iterator(world, sim_parameters, callback):
    '''
    Generic caller to call any function to be executed in each ecological 
    cell in sequence.
    
    @param world: dose_world.World object
    @param sim_parameters: simulation parameters dictionary (see Examples)
    @param callback: function to be executed
    @return: none
    '''
    for x in range(sim_parameters["world_x"]):
        for y in range(sim_parameters["world_y"]):
            for z in range(sim_parameters["world_z"]):
                if len(inspect.getargspec(callback)[0]) == 5:
                    callback(world, x, y, z)
                else:
                    callback(world)


def deploy_0(sim_parameters, populations, pop_name, world):
    '''
    Organism deployment scheme 0 - User defined deployment scheme. This is 
    called when "deployment_code" in simulation parameters = 0. User will 
    have to provide a deployment scheme/function as "deployment_scheme" in 
    simulation parameters.
    
    @param sim_parameters: simulation parameters dictionary (see Examples)
    @param populations: dictionary of population objects
    @param pop_name: population name
    @param world: dose_world.World object
    @return: none.
    '''
    sim_parameters["deployment_scheme"](populations, pop_name, world)


def deploy_1(sim_parameters, populations, pop_name, world):
    '''
    Organism deployment scheme 1 - Single ecological cell deployment 
    scheme. This is called when "deployment_code" in simulation parameters 
    = 1. In this scheme, all organisms in the specified population 
    (specified by population name) will be deployed in the ecological cell 
    specified in "population_locations" of simulation parameters.
    
    @param sim_parameters: simulation parameters dictionary (see Examples)
    @param populations: dictionary of population objects
    @param pop_name: population name
    @param world: dose_world.World object
    @return: none.
    '''
    position = sim_parameters["population_names"].index(pop_name)
    locations = [location for location in sim_parameters["population_locations"][position]]
    (x, y, z) = coordinates(locations[0])
    world.ecosystem[x][y][z]['organisms'] = sim_parameters["population_size"]
    for individual in populations[pop_name].agents:
        individual.status['location'] = locations[0]


def deploy_2(sim_parameters, populations, pop_name, world):
    '''
    Organism deployment scheme 2 - Random deployment scheme. This is called 
    when "deployment_code" in simulation parameters = 2. In this scheme, 
    all organisms in the specified population (specified by population 
    name) will be randomly deployed across the list of ecological cells 
    specified as "population_locations" of simulation parameters.
    
    @param sim_parameters: simulation parameters dictionary (see Examples)
    @param populations: dictionary of population objects
    @param pop_name: population name
    @param world: dose_world.World object
    @return: none
    '''
    position = sim_parameters["population_names"].index(pop_name)
    locations = [location for location in sim_parameters["population_locations"][position]]
    for individual in populations[pop_name].agents:
            location = random.choice(locations)
            (x, y, z) = coordinates(location)
            while world.ecosystem[x][y][z]['organisms'] >= sim_parameters["eco_cell_capacity"]:
                location = random.choice(locations)
                (x, y, z) = coordinates(location)
            world.ecosystem[x][y][z]['organisms'] = world.ecosystem[x][y][z]['organisms'] + 1
            individual.status['location'] = location


def deploy_3(sim_parameters, populations, pop_name, world):
    '''
    Organism deployment scheme 3 - Even deployment scheme. This is called 
    when "deployment_code" in simulation parameters = 3. In this scheme, 
    all organisms in the specified population (specified by population 
    name) will be evenly deployed (as much as possible) across the list of 
    ecological cells specified as "population_locations" of simulation 
    parameters.
    
    @param sim_parameters: simulation parameters dictionary (see Examples)
    @param populations: dictionary of population objects
    @param pop_name: population name
    @param world: dose_world.World object
    @return: none
    '''

    position = sim_parameters["population_names"].index(pop_name)
    locations = [location for location in sim_parameters["population_locations"][position]]
    iterator = 0
    for i in range(sim_parameters["population_size"]):
        individual = populations[pop_name].agents[i]
        location = locations[iterator]
        (x, y, z) = coordinates(location)
        world.ecosystem[x][y][z]['organisms'] += 1
        individual.status['location'] = location
        iterator += 1
        if iterator == len(locations):
            iterator = 0


def deploy_4(sim_parameters, populations, pop_name, world):
    '''
    Organism deployment scheme 4 - Centralized deployment scheme. This is 
    called when "deployment_code" in simulation parameters = 4. In this 
    scheme, all organisms in the specified population (specified by 
    population name) will be deployed onto the ecological cell specified 
    in "population_locations" of simulation parameters, up to the organism 
    capacity of the ecological cell (specified in "eco_cell_capacity" of 
    simulation parameters). In event that the specified deployment 
    locations is filled, undeployed organisms will be randomly deployed 
    onto adjacent cells.
    
    @param sim_parameters: simulation parameters dictionary (see Examples)
    @param populations: dictionary of population objects
    @param pop_name: population name
    @param world: dose_world.World object
    @return: none
    '''
    position = sim_parameters["population_names"].index(pop_name)
    locations = [location for location in sim_parameters["population_locations"][position]]
    adj_cells = adjacent_cells(sim_parameters, locations[0])
    for group in range((sim_parameters["population_size"] / sim_parameters["eco_cell_capacity"]) + 1):
        start = sim_parameters["eco_cell_capacity"] * group
        end = start + sim_parameters["eco_cell_capacity"]
        for x in range(start, end):
            if x == sim_parameters["population_size"]:
                break
            individual = populations[pop_name].agents[x]
            if x > (sim_parameters["eco_cell_capacity"] - 1):
                location = random.choice(adj_cells)
                (x, y, z) = coordinates(location)
                while world.ecosystem[x][y][z]['organisms'] > sim_parameters["eco_cell_capacity"]:
                    location = random.choice(adj_cells)
                    (x, y, z) = coordinates(random.choice(adj_cells))
            (x, y, z) = coordinates(location)
            world.ecosystem[x][y][z]['organisms'] += 1
            individual.status['location'] = location


def interpret_chromosome(sim_parameters, populations, pop_name, world):
    '''
    Function to call Ragaraja interpreter to express / execute the genome 
    for each organism in a population. The Turing tape (array) after 
    execution will be logged as "blood" of each organism. Ragaraja 
    interpreter will use the temporary_input and temporary_output lists 
    from each ecological cell as the input data and output data respectively 
    for genome execution - this will resemble consumption of environmental 
    resources and replenishing of environmental resources or dumping of 
    wastes respectively.
    
    @param sim_parameters: simulation parameters dictionary (see Examples)
    @param populations: dictionary of population objects
    @param pop_name: population name
    @param world: dose_world.World object
    @return: none
    '''

    chromosome = None

    for i in range(len(populations[pop_name].agents)):

        individual = populations[pop_name].agents[i]
        location = individual.status['location']
        (x, y, z) = coordinates(location)

        if sim_parameters["clean_cell"]:
            chromosome = [0] * sim_parameters["max_tape_length"]
        else:
            chromosome = populations[pop_name].agents[i].status['blood']
            if chromosome is None:
                chromosome = [0] * sim_parameters["max_tape_length"]

        for chromosome_count in range(len(individual.genome)):
            inputdata = world.ecosystem[x][y][z]['local_input']
            output = world.ecosystem[x][y][z]['local_output']
            source = ''.join(individual.genome[chromosome_count].sequence)
            chromosome = populations[pop_name].agents[i].status['blood']
            try:
                (chromosome, apointer, inputdata, output, source, spointer) = \
                    register_machine.interpret(source, ragaraja.ragaraja, 3, inputdata, chromosome,
                                               sim_parameters["max_tape_length"], sim_parameters["max_codon"])
            except Exception as e: 
                error_msg = '|'.join(['Error at Chromosome_' + str(chromosome_count), str(e)])
                populations[pop_name].agents[i].status['chromosome_error'] = error_msg
                populations[pop_name].agents[i].status['blood'] = chromosome

            populations[pop_name].agents[i].status['blood'] = chromosome
            world.ecosystem[x][y][z]['temporary_input'] = inputdata
            world.ecosystem[x][y][z]['temporary_output'] = output


def step(populations, pop_name, sim_functions):
    '''
    Performs a generational step for a population
        - Prepopulation control
        - Mutations
        - Before mating fitness measurement
        - Mating
        - Postpopulation control
        - Generational events
        - After mating fitness measurement
        - Generate a textual report for the current generation
    
    @param populations: dictionary of population objects
    @param pop_name: population name
    @param sim_functions: implemented simulation functions 
    (see dose.dose_functions)
    @return: report as a string
    '''
    if populations[pop_name].generation > 0:
        sim_functions.prepopulation_control(populations, pop_name)
    for organism in populations[pop_name].agents:
        sim_functions.mutation_scheme(organism)
    sim_functions.fitness(populations, pop_name)
    sim_functions.mating(populations, pop_name)
    sim_functions.postpopulation_control(populations, pop_name)
    sim_functions.generation_events(populations, pop_name)
    populations[pop_name].generation = populations[pop_name].generation + 1
    sim_functions.fitness(populations, pop_name)
    return sim_functions.population_report(populations, pop_name)


def report_generation(sim_parameters, populations, pop_name, sim_functions, generation_count):
    '''
    Performs a generational step (using step function) for a population 
    and writes out the resulting report into results text file.
    
    @param sim_parameters: simulation parameters dictionary (see Examples)
    @param populations: dictionary of population objects
    @param pop_name: population name
    @param sim_functions: implemented simulation functions (see 
    dose.dose_functions)
    @param generation_count: current generation count for reporting
    @return: none
    '''
    report = step(populations, pop_name, sim_functions)
    if generation_count % int(sim_parameters["fossilized_frequency"]) == 0:
        file = '%s%s_%s_' % (sim_parameters["directory"], sim_parameters["simulation_name"], pop_name)
        freeze_population(file, sim_parameters["fossilized_ratio"], populations, pop_name)
    if generation_count % int(sim_parameters["print_frequency"]) == 0:
        f = open(('%s%s_%s.result.txt' % (sim_parameters["directory"], sim_parameters["simulation_name"],
                                          pop_name)), 'a')
        dtstamp = str(datetime.utcnow())
        f.write('\n'.join(['\n' + dtstamp, 'GENERATION: ' + str(generation_count), str(report)]))
        f.write('\n')
        f.close


def bury_world(sim_parameters, world, generation_count):
    '''
    Function to bury entire world into a file.
    
    @param sim_parameters: simulation parameters dictionary (see Examples)
    @param world: dose_world.World object
    @param generation_count: current generation count for file name generation
    '''

    filename = '%s%s_gen%s.eco' % (sim_parameters["directory"], sim_parameters["simulation_name"],
                                   str(generation_count))
    f = open(filename, 'wb')
    pickle.dump(world, f)
    f.close()


def excavate_world(eco_file):
    '''
    Excavate buried world from file.
    
    @param eco_file: buried world file generated by bury_world function.
    @return: excavated dose_world.World object
    '''
    f = open(eco_file, 'rb')
    return pickle.load(f)


def freeze_population(file, proportion, populations, pop_name):
    '''
    Function to freeze part or whole of the population into a file. If 
    the number of organisms is less than 101, the entire population will 
    be frozen.
    
    @param file: file name prefix
    @param proportion: proportion of the population to freeze
    @param populations: dictionary of population objects
    @param pop_name: population name to freeze
    '''
    if proportion > 1.0: proportion = 1.0
    agents = populations[pop_name].agents
    if len(agents) < 101 or len(agents) * proportion < 101:
        sample = deepcopy(populations[pop_name])
    else:
        new_agents = [agents[random.randint(0, len(agents) - 1)] for x in range(int(len(agents) * proportion))]
        sample = deepcopy(populations[pop_name])
        sample.agents = new_agents
    name = ''.join([file, 'pop', str(populations[pop_name].generation), '_', str(len(sample.agents)), '.gap'])
    f = open(name, 'wb')
    pickle.dump(sample, f)
    f.close()


def revive_population(gap_file):
    '''
    Revive population from a frozen population file.
    
    @param gap_file: frozen population file generated by freeze_population
    function.
    @return: revived population object
    '''
    f = open(gap_file, 'rb')
    return pickle.load(f)


def write_parameters(sim_parameters, pop_name):
    '''
    Function to write simulation parameters into results text file as 
    header.
    
    @param sim_parameters: simulation parameters dictionary (see Examples)
    @param pop_name: population name
    @return: none
    '''
    f = open(('%s%s_%s.result.txt' % (sim_parameters["directory"], sim_parameters["simulation_name"], pop_name)), 'a')
    f.write("SIMULATION: %s \n" + ("-" * 70) % sim_parameters["simulation_name"])
    if ('database_source' in sim_parameters) or  ('sim_folder' in sim_parameters):
        f.write("SIMULATION REVIVAL STARTED: %s\n\n" % \
                sim_parameters["starting_time"])
    else:
        f.write("SIMULATION STARTED: %s\n\n" % sim_parameters["starting_time"])
    for key in sim_parameters:
        if key not in ('deployment_scheme', 'directory', 'sim_folder'):
            f.write("%s : %s\n" % (key, sim_parameters[key]))
    f.write("\n\nREPORT\n" + + ("-" * 70))


def close_results(sim_parameters, pop_name):
    '''
    Function to write footer of results text file.
    
    @param sim_parameters: simulation parameters dictionary (see Examples)
    @param pop_name: population name
    @return: none
    '''
    f = open(('%s%s_%s.result.txt' % (sim_parameters["directory"],
                                      sim_parameters["simulation_name"], 
                                      pop_name)), 'a')
    f.write("-" * 70 + "\nSIMULATION ENDED: " + str(datetime.utcnow()))
    f.close()
