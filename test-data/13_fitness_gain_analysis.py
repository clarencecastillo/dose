import sys, os, random
cwd = os.getcwd().split(os.sep)
cwd[-1] = 'dose'
cwd = os.sep.join(cwd)
sys.path.append(cwd)
import analytics

def analyze_fitness(fitness): return fitness

for trial in range(1, 26):
	analysis_name = "T" + str(trial) + "_FPS_11x0_gain4"
	sim13_analysis = analytics.Analysis(analysis_name + ".db", 
										"pop_01", starting_time = 'default')
	sim13_analysis.analyze_individual_status_by_generation(analysis_name + ".csv", 
														   "fitness", analyze_fitness, 
														   {"Average":analytics.average, 
															"STDEV":analytics.standard_deviation}, 
															'all')