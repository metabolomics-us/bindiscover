from subprocess import Popen
import os
import sys
from pprint import pprint

def divide_chunks(l, n):
     
    # looping till length l
    for i in range(0, len(l), n):
        yield l[i:i + n]

def build_one_command(temp_file_name):
    command_string=f'python3 ./code/generate_signifigance_test_matrices.py 0 1 {temp_file_name}\n'
    return command_string

if __name__=="__main__":

    min_fold_change=sys.argv[1]
    cores_available=int(sys.argv[2])


    os.system('mkdir -p ../results/'+str(min_fold_change)+'/step_6_b_generate_signifigance_test_matrices/')
    os.system('touch ../results/'+str(min_fold_change)+'/step_6_b_generate_signifigance_test_matrices/dummy.txt')


    full_file_list=os.listdir('../results/'+str(min_fold_change)+'/step_6_generate_fold_matrices/')
    full_file_list.remove('dummy.txt')

    full_file_list_list=list(divide_chunks(full_file_list,cores_available))

    full_command_list_list=list()
    for i in range(len(full_file_list_list)):
        full_command_list_list.append(list())
        for j in range(len(full_file_list_list[i])):
            full_command_list_list[i].append(
                build_one_command(full_file_list_list[i][j])
            )
    
    for i in range(len(full_command_list_list)):
        print('command set: '+str(i))
        temp_Popen_list=[Popen(element,shell=True) for element in full_command_list_list[i]]
        for command in temp_Popen_list:
            command.wait()