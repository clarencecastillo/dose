'''
Application Programming Interface (API) for DOSE (digital organism 
simulation environment). This contains the main functions and operations 
needed to write a DOSE simulation. This file will be imported as top 
level (from dose import *) when DOSE is imported; hence, all functions in 
this file can be assessed at top level.

Date created: 27th September 2013
'''
import sys, os, random, inspect

import ragaraja, register_machine
import dose_world
import genetic, simulation_calls

from simulation_calls import spawn_populations, simulation_core
from simulation_calls import excavate_world, revive_population

from database_calls import connect_database, db_reconstruct_simulation_parameters
from database_calls import db_reconstruct_population, db_reconstruct_world, db_list_simulations

class dose_functions():
    '''
    Abstract class to contain all of the simulation-specific functions 
    (functions that vary with each simulation) that are to be defined / 
    implemented by the user to be used in a simulation. This class should 
    be inherited by every simulation to over-ride each function / method. 
    
    This set of functions / methods is consolidated functions / methods to 
    be over-ridden from genetic.Organism, genetic.Population, and 
    dose_world.World classes. As a result, the functions / methods can be
    working at different levels - at the level of individual organisms, 
    at the level of entire population(s), or at the level of the world.
    
    Please see the examples in examples directory on its use.
    '''
    def mutation_scheme(self, organism):
        '''
        Method / function to trigger mutational events in each chromosome 
        of the genome within an organism. This function works at the 
        level of individual organisms.
        
        @param organism: genetic.Organism object
        @return: None
        '''
        raise NotImplementedError
    def prepopulation_control(self, Populations, pop_name):
        '''
        Method / function to trigger population control events before 
        mating event in each generation. For example, it can be used to
        simulate pre-puberty (childhood) death. This function works at 
        the level of entire population(s).
        
        @param Populations: A dictionary containing one or more populations 
        where the value is a genetic.Population object.
        @param pop_name: Name of the population which is used as key in 
        the the dictionary (Populations parameter).
        @return: None
        '''
        raise NotImplementedError
    def fitness(self, Populations, pop_name):
        '''
        Method / function to calculate the fitness score of each organism 
        within the population(s). This function works at the level of 
        entire population(s) even though fitness calculation occurs at the 
        organism level. The fitness of each organism may be stored in 
        Organism.status['fitness'] and may be used by mating scheme.
        
        @param Populations: A dictionary containing one or more populations 
        where the value is a genetic.Population object.
        @param pop_name: Name of the population which is used as key in 
        the the dictionary (Populations parameter).
        @return: None
        '''
        raise NotImplementedError
    def mating(self, Populations, pop_name):
        '''
        Method / function to trigger mating events in each generation. For 
        example, it can be used to simulates mate choices and progeny size. 
        A support function provided is genetic.crossover() function which 
        generates one random crossover operation between 2 chromosomes, to 
        simulate meiosis crossover. This function may also use one or more 
        of the dose.filter_XXX() functions to select or choose suitable 
        mates. This function works at the level of entire population(s), 
        which means that this function will have 
            - to manage mating scheme and progeny (offspring) generation 
            for the entire population
            - add or replace offsprings into the respective population(s)
            - (optional) store identity of parent(s) as list; for example, 
            [parentA identity, parentB identity], in offspring's 
            status['parents'] for ancestral tracing.
        
        @param Populations: A dictionary containing one or more populations 
        where the value is a genetic.Population object.
        @param pop_name: Name of the population which is used as key in 
        the the dictionary (Populations parameter).
        @return: None
        '''
        raise NotImplementedError
    def postpopulation_control(self, Populations, pop_name):
        '''
        Method / function to trigger population control events after 
        mating event in each generation. For example, it can be used to
        simulate old-age death. This function works at the level of entire 
        population(s).
        
        @param Populations: A dictionary containing one or more populations 
        where the value is a genetic.Population object.
        @param pop_name: Name of the population which is used as key in 
        the the dictionary (Populations parameter).
        @return: None
        '''
        raise NotImplementedError
    def generation_events(self, Populations, pop_name):
        '''
        Method / function to trigger other defined events in each 
        generation. For example, it can be used to simulate catastrophe 
        or epidemic that does not occur regularly, or simulates unusual 
        occurrences of multiple mutation events. This function works at 
        the level of entire population(s).
        
        @param Populations: A dictionary containing one or more populations 
        where the value is a genetic.Population object.
        @param pop_name: Name of the population which is used as key in 
        the the dictionary (Populations parameter).
        @return: None
        '''
        raise NotImplementedError
    def population_report(self, Populations, pop_name):
        '''
        Method / function to generate a text report of the population(s) 
        and/or each organisms within the population at regular intervals, 
        within the simulation, as determined by "print_frequency" in the 
        simulation parameters. This function works at the level of entire 
        population(s).
        
        @param Populations: A dictionary containing one or more populations 
        where the value is a genetic.Population object.
        @param pop_name: Name of the population which is used as key in 
        the the dictionary (Populations parameter).
        @return: Entire report of a population at a generation count in a 
        string. To make it human-readable, usually the report string is 
        both tab-delimited (for one organism) and newline-delimited (for 
        entire population). 
        '''
        raise NotImplementedError
    def organism_movement(self, Populations, pop_name, World):
        '''
        organism_movement and organism_location are both methods / 
        functions to execute movement of organisms within the world. The 
        semantic difference between organism_movement and organism_location 
        is that organism_movement is generally used for short travels 
        while organism_location is used for long travel. For example, 
        organism_movement can be used to simulate foraging or nomadic
        behaviour. This function works at both the level of entire 
        population(s) and world.
        
        For each organism to move, this function will have to 
            - update the number of organisms in each 
            World.ecosystem[x-axis][y-axis][z-axis]['organisms']
            - update the respective Organism's location in the status 
            dictionary (Population[pop_name].agents[<index>].status['location'])
        
        @param Populations: A dictionary containing one or more populations 
        where the value is a genetic.Population object.
        @param pop_name: Name of the population which is used as key in 
        the the dictionary (Populations parameter).
        @param World: dose_world.World object.
        @return: None
        '''
        raise NotImplementedError
    def organism_location(self, Populations, pop_name, World):
        '''
        organism_movement and organism_location are both methods / 
        functions to execute movement of organisms within the world. The 
        semantic difference between organism_movement and organism_location 
        is that organism_movement is generally used for short travels 
        while organism_location is used for long travel. For example, 
        organism_movement can be used to simulate long distance migration, 
        such as air travel. This function works at both the level of entire 
        population(s) and world.
        
        For each organism to move, this function will have to 
            - update the number of organisms in each 
            World.ecosystem[x-axis][y-axis][z-axis]['organisms']
            - update the respective Organism's location in the status 
            dictionary (Population[pop_name].agents[<index>].status['location'])
        
        @param Populations: A dictionary containing one or more populations 
        where the value is a genetic.Population object.
        @param pop_name: Name of the population which is used as key in 
        the the dictionary (Populations parameter).
        @param World: dose_world.World object.
        @return: None
        '''
        raise NotImplementedError
    def ecoregulate(self, World): 
        '''
        Method / function for broad spectrum management of the entire 
        ecosystem defined as World.ecosystem[x-axis][y-axis][z-axis]
        ['local_input'] and World.ecosystem[x-axis][y-axis][z-axis]
        ['local_output']). For example, it can be used to simulate 
        temperature, solar radiation, or resource gradients. This function 
        works at the level of the world.
        
        @param World: dose_world.World object.
        @return: None
        '''
        raise NotImplementedError
    def update_ecology(self, World, x, y, z):
        '''
        Method / function to process the input and output from the 
        activities of the organisms in the current ecological cell (defined 
        as World.ecosystem[x-axis][y-axis][z-axis]['temporary_input'] and 
        World.ecosystem[x-axis][y-axis][z-axis]['temporary_output']) into 
        a local ecological cell condition (defined as World.ecosystem
        [x-axis][y-axis][z-axis]['local_input'] and World.ecosystem[x-axis]
        [y-axis][z-axis]['local_output']), and update the
        ecosystem (which is essentially the ecological cell adjacent to 
        World.ecosystem[x-axis][y-axis][z-axis]). For example, it can be 
        used to simulate secretion of chemicals or use of resources (such 
        as food) by organisms, and diffusion of secretions to the 
        neighbouring ecological cells. 
        
        Essentially, this function simulates the "diffusion" of local 
        situation outwards. This function works at the level of the world.
        
        @param World: dose_world.World object.
        @param x: x-axis of the World.ecosystem to identify the cell.
        @param y: y-axis of the World.ecosystem to identify the cell.
        @param z: z-axis of the World.ecosystem to identify the cell.
        @return: None
        '''
        raise NotImplementedError
    def update_local(self, World, x, y, z):
        '''
        Method / function to update local ecological cell condition (defined 
        as World.ecosystem[x-axis][y-axis][z-axis]['local_input'] and 
        World.ecosystem[x-axis][y-axis][z-axis]['local_output']) from the 
        ecosystem (World.ecosystem [x-axis][y-axis][z-axis]['local_input'] 
        and World.ecosystem[x-axis][y-axis][z-axis]['local_output'] of 
        adjacent ecological cells). For example, it can be used to 
        simulates movement or diffusion of resources from the ecosystem to 
        local. 
        
        Essentially, this function is the reverse of update_ecology() 
        function. In this case, the local ecological cell is affected by 
        adjacent conditions. This function works at the level of the world.
        
        @param World: dose_world.World object.
        @param x: x-axis of the World.ecosystem to identify the cell.
        @param y: y-axis of the World.ecosystem to identify the cell.
        @param z: z-axis of the World.ecosystem to identify the cell.
        @return: None
        '''
        raise NotImplementedError
    def report(self, World):
        '''
        Method / function to generate a text report of the ecosystem status 
        (World.ecosystem) at regular intervals, within the simulation, as 
        determined by "print_frequency" in the simulation parameters. This 
        function works at the level of world.
        
        @param World: dose_world.World object.
        @return: Entire report of a population at a generation count in a 
        string. To make it human-readable, usually the report string is 
        both tab-delimited (for one organism) and newline-delimited (for 
        entire population). 
        '''
        raise NotImplementedError
    def deployment_scheme(self, Populations, pop_name, World):
        '''
        Method / function to implement a user-specific / simulation-
        specific deployment scheme used to deploy organisms into the 
        World. This function will only be used when "deployment_code" in 
        simulation parameters dictionary equals to 0. This function works 
        at all three levels - organism(s), population(s), and world.
        
        @param Populations: A dictionary containing one or more populations 
        where the value is a genetic.Population object.
        @param pop_name: Name of the population which is used as key in 
        the the dictionary (Populations parameter).
        @param World: dose_world.World object.
        @return: None
        '''
        raise NotImplementedError
    def database_report(self, con, cur, start_time,
                        Population, World, generation_count):
        '''
        Method / function to implement database logging of each organism 
        in each population, and the ecosystem status. The frequency of 
        logging is determined by "database_logging_frequency" in the 
        simulation_parameters.
        
        Three database tables had been defined for use in this function:
            - organisms where the structure is (start_time text, pop_name 
            text, org_name text, generation text, key text, value text)
            - world where the structure is (start_time text, x text, 
            y text, z text, generation text, key text, value text)
            - miscellaneous where the structure is (start_time text, 
            generation text, key text, value text)
            
        The logical purpose of these tables are:
            - organisms, to log status of each organism. Starting time of 
            current simulation (start_time), population name (pop_name), 
            organism name (org_name), and generation count (generation) 
            are used as complex primary key to identify the organism 
            within a specific simulation at a specific generation. Key and 
            value pair makes up the actual data to be logged where key is 
            the field and value is the datum or attribute.
            - world, to log the status of the ecosystem. Starting time of 
            current simulation (start_time), location of ecological cell 
            (x, y, z as coordinates), and generation count (generation) 
            are used as complex primary key to identify the ecological cell
            within a specific simulation at a specific generation. Key and 
            value pair makes up the actual data to be logged where key is 
            the field and value is the datum or attribute.
            - miscellaneous, is used to log other undefined data. Starting 
            time of current simulation (start_time), and generation count 
            (generation) are used as complex primary key to identify a 
            specific generation within a specific simulation. Key and 
            value pair makes up the actual data to be logged where key is 
            the field and value is the datum or attribute.
        
        @param con: Database connector. See Python DB-API for details.
        @param cur: Database cursor. See Python DB-API for details.
        @param start_time: Starting time of current simulation in the 
        format of <date>-<seconds since epoch>; for example, 
        2013-10-11-1381480985.77.
        @param Population: A dictionary containing one or more populations 
        where the value is a genetic.Population object.
        @param World: dose_world.World object.
        @param generation_count: Current number of generations simulated.
        @return: None
        '''
        raise NotImplementedError

class event_broker():

    def __init__(self):
        self.observers = {}

    def log(self, message):
        print(message)

    def publish(self, event, data):
        if event in self.observers:
            for callback in self.observers[event]:
                callback(data)

    def subscribe(self, event, callback):
        if event in self.observers:
            self.observers[event].append(callback)
        else:
            self.observers[event] = [callback]

def database_report_populations(con, cur, start_time, 
                                Populations, generation_count):
    '''
    Function to log organisms' status and genome into database. Organisms' 
    status is implemented as a dictionary and each key-value pair in the 
    status is logged as a separate record. Similarly, each chromosome is 
    logged as a separately record. A combination of starting time of the 
    simulation, population name, organism's name, and the generation can 
    identify all the data specific to an organism within a population, at 
    a specific generation within a simulation. This function is a complete 
    logger - it logs everything there is about an organism. There is 
    nothing to log for populations. 
    
    The following transformations of data are made:
        - Organism.status['blood'] is a list of numbers resulting from 
        interpreting the genome by Ragaraja interpreter. The numbers are 
        concatenated and delimited by '|'. For example, [1, 2, 3] ==> 1|2|3
        - Organism.status['location'] is a tuple of 3 integers (x, y, z) for
        location of ecological cell. The numbers are concatenated and 
        delimited by '|'. For example, (2, 3, 4) ==> 2|3|4
        - Each chromosome in Organism.genome is a list of bases. These bases 
        are concatenated with no delimiter. For example, [1, 2, 3] ==> 123
    
    @param con: Database connector. See Python DB-API for details.
    @param cur: Database cursor. See Python DB-API for details.
    @param start_time: Starting time of current simulation in the 
    format of <date>-<seconds since epoch>; for example, 
    2013-10-11-1381480985.77.
    @param Populations: A dictionary containing one or more populations 
    where the value is a genetic.Population object.
    @param generation_count: Current number of generations simulated.
    @return: None
    '''
    generation = str(generation_count)
    for pop_name in list(Populations.keys()):
        for org in Populations[pop_name].agents:
            org_name = str(org.status['identity'])
            # log each item in Organism.status dictionary
            for key in [key for key in list(org.status.keys())
                               if key != 'identity']:
                if key in ('blood', 'location', 'parents'):
                    try:
                        # TypeError will occur when genome/chromosomes are 
                        # not interpreted
                        value = '|'.join([str(x) for x in org.status[key]])
                    except TypeError: value = ''
                else: 
                    value = str(org.status[key])
                cur.execute('insert into organisms values (?,?,?,?,?,?)', 
                    (str(start_time), str(pop_name), org_name, 
                     generation, key, value))
            # log each chromosome sequence
            for chromosome_count in range(len(org.genome)):
                key = 'chromosome_' + str(chromosome_count)
                sequence = ''.join(org.genome[chromosome_count].sequence)
                cur.execute('insert into organisms values (?,?,?,?,?,?)', 
                    (str(start_time), str(pop_name), org_name, 
                     generation, key, sequence))
    con.commit()
    
def database_report_world(con, cur, start_time, World, generation_count):
    '''
    Function to log World.ecosystem into database. The ecosystem is made 
    up of a collection of ecological cells, identified by the (x, y, z) 
    coordinates within the ecosystem. Each ecological cell is implemented 
    as a dictionary of ecological status. This function logs the entire 
    set of ecological cells; thus, a complete logger. Each ecosystem within 
    a specific simulation, at a specific time/generation, can be identified 
    by a combination of starting time of simulation and generation. Each 
    ecological cell within the ecosystem can be further identified by the 
    (x, y, z) coordinates.
    
    @param con: Database connector. See Python DB-API for details.
    @param cur: Database cursor. See Python DB-API for details.
    @param start_time: Starting time of current simulation in the 
    format of <date>-<seconds since epoch>; for example, 
    2013-10-11-1381480985.77.
    @param World: dose_world.World object.
    @param generation_count: Current number of generations simulated.
    @return: None
    '''
    generation = str(generation_count)
    ecosystem = World.ecosystem
    location = [(x, y, z) 
                for x in range(len(ecosystem))
                    for y in range(len(ecosystem[x]))
                        for z in range(len(ecosystem[x][y]))]
    for cell in location:
        eco_cell = ecosystem[cell[0]][cell[1]][cell[2]]
        for key in list(eco_cell.keys()):
            value = str(eco_cell[key])
            cur.execute('insert into world values (?,?,?,?,?,?,?)', 
                        (str(start_time), 
                         str(cell[0]), str(cell[1]), str(cell[2]), 
                         generation, key, value))
    con.commit()
    
def filter_deme(deme_name, agents):
    '''
    Function to identify organisms (agents) with a specific sub-population 
    name (also known as deme) within a population. Demes can be considered 
    as strains in bacteria, breed or sub-species in animals, races in 
    humans, and cultivars in plants. This function is can be used to 
    support the identification of suitable mates for mating schemes.
    
    @param deme_name: Name of deme (sub-population name)
    @type deme_name: string
    @param agents: A list of organisms, such as Population.agents.
    @return: List of Organism objects
    '''
    extract = [individual for individual in agents
               if individual.status['deme'].upper() == deme_name.upper()]
    return extract
    
def filter_gender(gender, agents):
    '''
    Function to identify organisms (agents) with a specific gender within 
    a population. This function is can be used to support the identification 
    of suitable mates for mating schemes.
    
    @param gender: Gender 
    @type gender: string
    @param agents: A list of organisms, such as Population.agents.
    @return: List of Organism objects
    '''
    extract = [individual for individual in agents
               if individual.status['gender'].upper() == gender.upper()]
    return extract

def filter_age(minimum, maximum, agents):
    '''
    Function to identify organisms (agents) within a certain age range in 
    a population. This function is can be used to support the identification 
    of suitable mates for mating schemes.
    
    @param minimum: Minimum age
    @type minimum: float
    @param maximum: Maximum age
    @type maximum: float
    @param agents: A list of organisms, such as Population.agents.
    @return: List of Organism objects
    '''
    extract = [individual for individual in agents
               if float(individual.status['age']) > (float(minimum) - 0.01) \
               and float(individual.status['age']) < float(maximum) + 0.01]
    return extract

def filter_location(location, agents):
    '''
    Function to identify organisms (agents) of a population within a 
    specific ecological cell. This function is can be used to support the 
    identification of suitable mates for mating schemes.
    
    @param location: (x, y, z) coordinates within the World.ecosystem
    @param location: tuple
    @param agents: A list of organisms, such as Population.agents.
    @return: List of Organism objects
    '''
    extract = [individual for individual in agents
               if individual.status['location'] == location]
    return extract

def filter_vitality(minimum, maximum, agents):
    '''
    Function to identify organisms (agents) within a certain vitality score 
    in a population. This function is can be used to support the identification 
    of suitable mates for mating schemes.
    
    @param minimum: Minimum vitality score
    @type minimum: float
    @param maximum: Maximum vitality score
    @type maximum: float
    @param agents: A list of organisms, such as Population.agents.
    @return: List of Organism objects
    '''
    extract = [individual for individual in agents
            if float(individual.status['vitality']) > (float(minimum) - 0.01) \
            and float(individual.status['vitality']) < float(maximum) + 0.01]
    return extract

def filter_status(status_key, condition, agents):
    '''
    Generic function to identity organisms (agents) within a population 
    via the status of organisms (Organism.status dictionary). This function 
    is can be used to support the identification of suitable mates for 
    mating schemes.
    
    @param status_key: Status (key in Organism.status dictionary) to filter
    @type status_key: string
    @param condition: Condition of the status (value in Organism.status 
    dictionary) to filter. This can be a unique condition, such as "True", 
    or a range, such as (minimum, maximum) to define the minimum and 
    maximum value of the condition.
    @param agents: A list of organisms, such as Population.agents.
    @return: List of Organism objects
    '''
    if type(condition) in (str, int, float, bool):
        extract = [individual for individual in agents 
                   if individual.status[status_key] == condition]
    else: 
        extract = [individual for individual in agents
        if float(individual.status[status_key]) > float(condition[0]) - 0.01 \
        and float(individual.status[status_key]) < float(condition[1]) + 0.01]
    return extract

def revive_simulation(rev_parameters, sim_functions, eb = event_broker()):

    eb.log("[" + rev_parameters["simulation_name"].upper() + " REVIVAL SIMULATION]")

    populations = {}
    world = None
    eb.log("Accessing simulation files directory")
    if "sim_folder" in rev_parameters:

        eb.log('Excavating World entity: ' + rev_parameters['eco_file'])
        world = excavate_world(rev_parameters['sim_folder'] + rev_parameters['eco_file'])

        eb.log('Updating parameters with World dimensions')
        rev_parameters["world_z"] = len(world.ecosystem[0][0][0])
        rev_parameters["world_y"] = len(world.ecosystem[0][0])
        rev_parameters["world_x"] = len(world.ecosystem[0])

        for i in range(len(rev_parameters["pop_files"])):
            eb.log('Reviving population file: ' + rev_parameters["pop_files"][i])
            pop_file = rev_parameters["sim_folder"] + rev_parameters["pop_files"][i]
            populations[rev_parameters["population_names"][i]] = revive_population(pop_file)

        eb.log('Updating revival generation start in simulation parameters')
        rev_parameters["rev_start"] = [populations[pop_name].generation for pop_name in populations]

    elif "database_source" in rev_parameters:

        eb.log('Constructing database directory')
        dbpath = os.sep.join([os.getcwd(), 'Simulations', rev_parameters["database_source"]])

        eb.log('Connecting to database file: ' + rev_parameters["database_source"])
        (con, cur) = connect_database(dbpath, None)

        if rev_parameters["simulation_time"] == 'default':
            eb.log('Acquiring simulation starting time')
            rev_parameters["simulation_time"] = db_list_simulations(cur)[0][0]

        eb.log('Reconstructing old simulation parameters')
        temp_parameters = db_reconstruct_simulation_parameters(cur, rev_parameters["simulation_time"])

        eb.log('Assimilating old simulation parameters with new simulation parameters')
        for key in temp_parameters:
            if key not in rev_parameters:
                rev_parameters[key] = temp_parameters[key]

        eb.log('Reconstructing World entity')
        world = db_reconstruct_world(cur, rev_parameters["simulation_time"], rev_parameters["rev_start"][0])

        eb.log('Updating population names parameter')
        for pop_name in rev_parameters["population_names"]:
            eb.log('Reconstructing population: ' + pop_name)
            populations[pop_name] = db_reconstruct_population(cur, rev_parameters["simulation_time"], pop_name,
                                                              rev_parameters["rev_start"][rev_parameters
                                                              ["population_names"].index(pop_name)])
        eb.log('Terminating database connection')
        con.close()

    eb.log('Updating last generation revival and population size simulation parameters')
    rev_parameters["rev_finish"] = [(populations[pop_name].generation + rev_parameters["extend_gen"])
                                    for pop_name in populations]
    rev_parameters["rev_pop_size"] = [len(populations[pop_name].agents) for pop_name in populations]

    eb.log('Starting simulation core')
    simulation_core(sim_functions, rev_parameters, populations, world)

def simulate(sim_parameters, sim_functions, eb = event_broker()):
    '''
    Function called by simulation to run the actual simulation based on a 
    set of parameters and functions.
    
    Simulation parameters dictionary contains the following keys:
        - simulation_name: Short name of the simulation.
        - population_names: Population name(s) in a list.
        - population_locations: A list of lists containing (x, y, z) 
        coordinates in the world to deploy the populations. There must be
        equal numbers of items in this list as population_names. However,
        one population can be deployed into one or more ecological cell(s).
        - deployment_code: Defines the type of deployment. Allowable values 
        are 0 (custom deployment scheme which users can implement by over-
        riding dose.dose_functions.deployment_scheme() function), 1 (all 
        organisms are in one location), 2 (organisms are randomly deployed 
        across a list of population locations given for the population), 3 
        (oganisms are randomly deployed across a list of population locations 
        given for the population), 4 (organisms are dispersed from a specific 
        location where the defined location will have the maximum allocation 
        before dispersal happens). (see simulation_calls.deployment for more 
        details)
        - chromosome_bases: List containing allowable bases in the 
        chromosome.
        - background_mutation: Defines background mutation rate where 0.01 
        represents 1% mutation and 0.5 represents 50% mutation.
        - additional_mutation: Defines mutation rate on top of background 
        mutation rate where 0.01 represents 1% mutation and 0.5 represents 
        50% mutation.
        - mutation_type: Type of mutation for background mutation rate. 
        Allowable values are 'point' (point mutation), 'insert' (insert a 
        base), 'delete' (delete a base), 'invert' (invert a stretch of the 
        chromosome), 'duplicate' (duplicate a stretch of the chromosome), 
        'translocate' (translocate  a stretch of chromosome to another 
        random position). 
        - chromosome_size: Number of bases for a chromosome.
        - genome_size: Number of chromosome(s) in a genome.
        - max_tape_length: Maximum number of cells for array used in 
        Ragaraja interpreter, which will be stored as 
        Organism.status['blood'].
        - clean_cell: Toggle to define if a new tape will be used for 
        Ragaraja interpreter. If True, at every Ragaraja interpretation of 
        the source (genome), a new tape ([0] * max_tape_length) will be 
        provided. If False, the tape stored as Organism.status['blood'] 
        (results from previous interpretation of genome) will be used.
        - interpret_chromosome: Toggle to define whether genomes will be 
        interpreted by Ragaraja. If True, chromosomes of each organism will 
        be interpreted by Ragaraja. 
        - max_codon: Maximum number of Ragaraja instructions to execute per 
        organism, which can be used as force-break from endless loop.
        - population_size: Number of organisms per population for initial 
        deployment.
        - eco_cell_capacity: Maximum number of organisms per ecological 
        cell at the time of deployment prior to the start of simulation.
        - world_x: Number of ecological cells in the x-axis of the 
        World.ecosystem.
        - world_y: Number of ecological cells in the y-axis of the 
        World.ecosystem.
        - world_z: Number of ecological cells in the z-axis of the 
        World.ecosystem.
        - goal: Goal for population to reach. This provides a goal for use 
        in fitness functions.
        - maximum_generations: Number of generations to simulate.
        - fossilized_ratio: Proportion (less than 1) of the population to 
        fossilize or freeze.
        - fossilized_frequency: Number of generations (intervals) for each 
        population fossilization or freezing event.
        - print_frequency: Number of generations (intervals) for each 
        reporting event into files.
        - ragaraja_version: Ragaraga instruction version to activate. (see 
        ragaraja.activate_version() for more details)
        - ragaraja_instructions: A list defining a set of Ragaraga 
        instructions to activate. This is only useful when ragaraja_version 
        is 0.
        - eco_buried_frequency: Number of generations (intervals) for each 
        ecological freezing or burial event.
        - database_file: Logging database file name. This file will be 
        located in <current working directory>/Simulations folder
        - database_logging_frequency: Number of generations (intervals) for 
        each database logging event. 
    
    Methods / Functions from dose.dose_functions class to be over-ridden 
    as simulation_functions (for more details, please look at 
    dose.dose_functions class):
        - organism_movement: organism_movement and organism_location are 
        both methods / functions to execute movement of organisms within 
        the world. The semantic difference between organism_movement and 
        organism_location is that organism_movement is generally used for 
        short travels while organism_location is used for long travel.
        - organism_location: organism_movement and organism_location are 
        both methods / functions to execute movement of organisms within 
        the world. The semantic difference between organism_movement and 
        organism_location is that organism_movement is generally used for 
        short travels while organism_location is used for long travel.
        - ecoregulate: Broad spectrum management of the entire 
        ecosystem.
        - update_ecology: Process the input and output from the activities 
        of the organisms in the current ecological cell into a local 
        ecological cell condition, and update the ecosystem.
        - update_local: Update local ecological cell condition from the 
        ecosystem.
        - report: Generate a text report of ecosystem (World.ecosystem) at 
        regular intervals, within the simulation, as determined by 
        "print_frequency" in the simulation parameters.
        - fitness: Calculate the fitness score of each organism within the 
        population(s).
        - mutation_scheme: Trigger mutational events in each chromosome 
        of the genome within an organism.
        - prepopulation_control: Trigger population control events before 
        mating event in each generation.
        - mating: Trigger mating events in each generation.
        - postpopulation_control: Trigger population control events after 
        mating event in each generation.
        - generation_events: Trigger other defined events in each 
        generation.
        - population_report: Generate a text report of the population(s) 
        and/or each organisms within the population at regular intervals, 
        within the simulation, as determined by "print_frequency" in the 
        simulation parameters
        - database_report: Implement database logging of each organism 
        in each population, and the ecosystem status. The frequency of 
        logging is determined by "database_logging_frequency" in the 
        simulation_parameters.
        - deployment_scheme: Implement a user-specific / simulation-
        specific deployment scheme used to deploy organisms into the 
        World. This function will only be used when "deployment_code" in 
        simulation parameters dictionary equals to 0.
    
    @param sim_parameters: Dictionary of simulation parameters
    @param sim_functions: A class inherited from dose.dose_functions
    class to implement all the needed simulation functions.
    @param eb: An instance of dose.event_broker
    '''

    eb.log("[" + sim_parameters["simulation_name"].upper() + " SIMULATION]")

    if "initial_chromosome" not in sim_parameters:
        eb.log('Adding initial chromosome to simulation parameters')
        sim_parameters["initial_chromosome"] = ['0'] * sim_parameters["chromosome_size"]

    eb.log('Adding deployment scheme to simulation parameters')
    sim_parameters["deployment_scheme"] = sim_functions.deployment_scheme

    eb.log('Constructing World entity')
    world = dose_world.World(sim_parameters["world_x"], sim_parameters["world_y"], sim_parameters["world_z"])
    eb.log('Spawning populations')
    populations = spawn_populations(sim_parameters)
    eb.log('Starting simulation core')
    simulation_core(sim_functions, sim_parameters, eb, populations, world)
