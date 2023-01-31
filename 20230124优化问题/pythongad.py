# %%
import pygad
import numpy as np
import matlab
import matlab.engine


# %%
#%% 2. 定义函数
eng = matlab.engine.start_matlab()
def fitness_func(solution, solution_idx):
# def black_box_function(mTMD1,mTMD2,fTMD1,fTMD2,fTMD3,dTMD1,dTMD2,dTMD3,xTMD1,xTMD2,xTMD3):
    fTMD1=solution[0]
    dTMD1=solution[1]
    xTMD1=solution[2]
    total_tmd_mass_ratio = 0.02 # 总质量比 The total mass ratio
    mass_six_span = 10007779.7 # 深中通道非通航桥六跨连续梁质量 The mass of 6-span continuous beam of the non-navigational bridge of the Zhenzhong-Link
    total_tmd_mass = total_tmd_mass_ratio * mass_six_span # 总质量 The total mass
    
    mTMD1=total_tmd_mass / number_of_tmds/6 # 质量 The mass mTMD1
    fTMD1=0.7+fTMD1*0.3 # 频率 The frequency fTMD1
    dTMD1=0.05+dTMD1*0.15 # 阻尼比 The damping ratio dTMD1
    xTMD1=xTMD1*660 # TMD1的x坐标 The x-coordinate of TMD1
    t_length=matlab.double(100) # 时间长度 The time length
    number_of_modes_to_control=matlab.double([1]) # 控制模态 The controlled modes
    # number_of_modes_to_control=matlab.double([1]) # 控制模态 The controlled modes
    number_of_modes_to_consider=10 # 考虑模态 The considered modes
    number_of_tmds=1 # TMD数量 The number of TMDs
    modal_damping_ratios=np.ones((1,number_of_modes_to_consider))*0.003 # 模态阻尼比 The modal damping ratios
    
    fitness = -eng.c_mode1(number_of_modes_to_control,number_of_modes_to_consider,number_of_tmds,modal_damping_ratios,t_length,mTMD1,fTMD1,dTMD1,xTMD1)
    # fitness = 1
    return fitness

# %%
gene_space = [{'low': 0, 'high': 1}]

fitness_function = fitness_func

num_generations = 1000 # Number of generations.
num_parents_mating = 7 # Number of solutions to be selected as parents in the mating pool.

# To prepare the initial population, there are 2 ways:
# 1) Prepare it yourself and pass it to the initial_population parameter. This way is useful when the user wants to start the genetic algorithm with a custom initial population.
# 2) Assign valid integer values to the sol_per_pop and num_genes parameters. If the initial_population parameter exists, then the sol_per_pop and num_genes parameters are useless.

sol_per_pop = 50 # Number of solutions in the population.
num_genes = 3 # Number of genes in the solution.


# num_generations = 10 # Number of generations.
# num_parents_mating = 7 # Number of solutions to be selected as parents in the mating pool.
# sol_per_pop = 11 # Number of solutions in the population.
# num_genes = 11 # Number of genes in the solution.

gene_space = [{'low': 0, 'high': 1}] * num_genes

last_fitness = 0
def callback_generation(ga_instance):
    global last_fitness
    print("Generation = {generation}".format(generation=ga_instance.generations_completed))
    print("Fitness    = {fitness}".format(fitness=ga_instance.best_solution()[1]))
    print("Change     = {change}".format(change=ga_instance.best_solution()[1] - last_fitness))
    last_fitness = ga_instance.best_solution()[1]

# Creating an instance of the GA class inside the ga module. Some parameters are initialized within the constructor.
ga_instance = pygad.GA(num_generations=num_generations,
                       num_parents_mating=num_parents_mating, 
                       fitness_func=fitness_function,
                       sol_per_pop=sol_per_pop, 
                       num_genes=num_genes,
                       gene_space=gene_space,
                       on_generation=callback_generation,
                       save_solutions=True)

print("Initial Population")
print(ga_instance.initial_population)

# Running the GA to optimize the parameters of the function.
ga_instance.run()

# # After the generations complete, some plots are showed that summarize the how the outputs/fitenss values evolve over generations.
# ga_instance.plot_fitness()

# # Returning the details of the best solution.
# solution, solution_fitness, solution_idx = ga_instance.best_solution()
# print("Parameters of the best solution : {solution}".format(solution=solution))
# print("Fitness value of the best solution = {solution_fitness}".format(solution_fitness=solution_fitness))
# print("Index of the best solution : {solution_idx}".format(solution_idx=solution_idx))



# if ga_instance.best_solution_generation != -1:
#     print("Best fitness value reached after {best_solution_generation} generations.".format(best_solution_generation=ga_instance.best_solution_generation))

# # Saving the GA instance.
# filename = 'genetic' # The filename to which the instance is saved. The name is without extension.
# ga_instance.save(filename=filename)

# # Loading the saved GA instance.
# loaded_ga_instance = pygad.load(filename=filename)
# loaded_ga_instance.plot_fitness()


