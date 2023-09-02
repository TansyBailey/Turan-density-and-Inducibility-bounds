###########
import math
###########

# Choose directory to save output text files to (uncomment below):
# savedirectory = ###


#####################################################


# Function to generate all unlabelled F-free graphs
def F_free_rgraphs(r,k,F):

    G = []

    for edge_size in range(math.comb(k,r)+1):
        G = G + list(hypergraphs.nauty(edge_size, k, uniform=r))

    G_structures = [IncidenceStructure(k,g) for g in G]
    
    F_free = [g for g in G_structures 
              if sum([int(sum(1 for _ in g.isomorphic_substructures_iterator(f)) > 0) 
              for f in F]) == 0]

    return(F_free)


######################################################


# Function to determine labelled isomorphism
def flag_isomorphic(F,G):

    # Check number of vertices:
    if F.num_points() != G.num_points():
        return False
     
    # Check type labels are identical:
    F_type_labels = [v for v in F.ground_set() if isinstance(v,str)]
    G_type_labels = [v for v in G.ground_set() if isinstance(v,str)]
    if sorted(F_type_labels) != sorted(G_type_labels):
        return False

    # Check if there is an isomorphism that preserves the labelling of the type vertices:
    type_dict = {v:v for v in F_type_labels}
    for iso in F.isomorphic_substructures_iterator(G,induced=True):
        if type_dict.items() <= iso.items(): 
            return True

    # If no isomorphism:  
    return(False)


######################################################


# Function to find probability p(F,G;H)
def flag_density(F,G,H):

    # Check F,G have the same labels:
    F_labels = [lab for lab in F.ground_set() if isinstance(lab, str)]
    G_labels = [lab for lab in G.ground_set() if isinstance(lab, str)]
    if F_labels.sort() != G_labels.sort():
        return 0

    # Check F,G have the same type:
    F_type = F.induced_substructure(F_labels)
    G_type = G.induced_substructure(G_labels)
    if not(flag_isomorphic(F_type,G_type)):
        return 0
    
    # Check graph dimensions make sense:
    type_len = len(F_labels)
    if (H.num_points() - type_len == 0 and 
        (F.num_points() - type_len > 0 
         or G.num_points() - type_len > 0)):
        return 0

    # Check H induces a structure isomorphic to the type:
    type_lst = set([tuple(sorted(list(d.values()))) 
                    for d in H.isomorphic_substructures_iterator(F_type,True)])
    if len(type_lst) == 0:
        return 0
    induced_types = [list(x) for x in type_lst]

    # Determine flags with underlying structure H:
    flag_probabilities = []
    num_F_unlabelled = F.num_points() - type_len 
    num_G_unlabelled = G.num_points() - type_len
    for label_lst in induced_types:
        label_perms = Permutations(label_lst,type_len)
        unlabelled = [v for v in H.ground_set() if v not in label_lst]
        unlabelled_dict = {v:v for v in unlabelled}
        dicts = [{p[i]:F_labels[i] for i in range(type_len)} for p in label_perms]
        for d in dicts: 
            H_copy = IncidenceStructure(H.ground_set(),H.blocks())
            H_copy.relabel(d|unlabelled_dict)
            unlabelled_perms = Permutations(unlabelled,num_F_unlabelled + num_G_unlabelled)
            perms_split = [[sorted(p[:num_F_unlabelled]),sorted(p[num_F_unlabelled:])] 
                           for p in unlabelled_perms]
            done_combinations = []
            num_induced_FG = 0 
            for lsts in perms_split:
                if lsts not in done_combinations:
                    done_combinations.append(lsts)
                    if (flag_isomorphic(H_copy.induced_substructure(F_labels + lsts[0]),F) and 
                        flag_isomorphic(H_copy.induced_substructure(F_labels + lsts[1]),G)):
                        num_induced_FG += 1
            flag_probabilities.extend(num_induced_FG/len(done_combinations))
    
    return(sum(flag_probabilities)/(len(induced_types)*math.factorial(type_len)))


######################################################


# Main function to find bound on desired quantity
def main_function(r,k,F,P,problem):

    # Generate unlabelled graphs of order k:
    k_graphs_list = F_free_rgraphs(r,k,F)
    num_graphs = len(k_graphs_list)

    # Save graph outputs to files:
    with open(savedirectory + "F-free graphs.txt", "w") as f:
        print('The edge sets of the F-free r-graphs on k vertices are:',file = f)
        for graph in k_graphs_list:
            print(graph.blocks(), file = f)
    with open(savedirectory + "output.txt", "w") as f:
        if num_graphs == 1:
            print('There is ' + 1 + " unlabelled " + str(r) + "-graph of order " + str(k) + '.',file = f)
        else:
            print('There are ' + str(num_graphs) + " unlabelled " + str(r) + "-graphs of order " + str(k) + '.',file = f)


    # Generate all types and flags:
    orders = [n for n in range(1, k-1) if n % 2 == k % 2]
    all_flags = []

    for v in orders:

        graphs_order_v = F_free_rgraphs(r,v,F)

        label_list = ['v' + str(i) for i in range(1, v+1)]
        [IS.relabel(label_list) for IS in graphs_order_v]

        flags_order_v = [[] for x in range(len(graphs_order_v))]
        all_flags_order_v = []
        num_t_flag_verts = int((k + v)/2)
        unlabelled_flags = F_free_rgraphs(r,num_t_flag_verts,F)

        vertices = list(range(0,num_t_flag_verts))
        type_permutations = Permutations(vertices,v)

        for perm in type_permutations:
            relabel_dict = {perm[i]:label_list[i] for i in range(v)} | {j:j for j in vertices if j not in perm}
            for graph in unlabelled_flags:
                flag = IncidenceStructure(graph.ground_set(),graph.blocks())
                flag.relabel(relabel_dict)
                which_type = [flag_isomorphic(flag.induced_substructure(label_list),g) for g in graphs_order_v]
                if any(which_type) and not any([flag_isomorphic(flag,x) for x in all_flags_order_v]):
                    flags_order_v[which_type.index(True)].append(flag)
                    all_flags_order_v.append(flag)
        all_flags.extend(all_flags_order_v)

        # Save type outputs to file:
        with open(savedirectory + "output.txt", "a") as f:
            if len(graphs_order_v) == 1 and v == 1:
                flag_lengths = [len(x) for x in flags_order_v]
                print('There is ' + str(len(graphs_order_v)) + ' type of order ' + str(v) 
                      + ', and ' + str(flag_lengths) + ' flags of order ' + str(num_t_flag_verts) + '.' ,file = f) 
            else:
                flag_lengths = [len(x) for x in flags_order_v]
                print('There are ' + str(len(graphs_order_v)) + ' types of order ' + str(v) 
                      + ', and ' + str(flag_lengths) + ' flags of order ' + str(num_t_flag_verts) + '.' ,file = f)

    # Save flag outputs to file:
    with open(savedirectory + "flags.txt", "w") as f:
        print('\nThe vertex and edge sets of the flags are:\n',file = f)
        for flag in all_flags:
            print('Vertex set: ', flag.ground_set(),';', 'Edge set: ', flag.blocks(), file = f)


    # Generate probability matrices:
    mat_dim = len(all_flags)
    H_matrices = [None] * num_graphs

    for H in k_graphs_list:
        H_matrix = matrix(RR,mat_dim,mat_dim)
        for F in all_flags:
            for G in all_flags:
                H_matrix[all_flags.index(F), all_flags.index(G)] = flag_density(F,G,IncidenceStructure(H.ground_set(),H.blocks()))
        H_matrices[k_graphs_list.index(H)] = H_matrix

        # Save probability matrices to file:
        with open(savedirectory + "probability matrices.txt", "w") as f:
            print('\nThe probability matrices for each unlabelled graph are:\n',file = f)
            for Hmatrix in H_matrices:
                print(Hmatrix,file = f)
                print('',file = f)


    # Inducibility problem
    if problem == 'I':

        # Written output:
        c_values = []
        for H_i in k_graphs_list:
            c_i = 0
            for P_i in P:
                p_sum = 0
                num_p_verts = len(P_i.ground_set())
                H_combos = Combinations(H_i.ground_set(),num_p_verts)
                for comb in H_combos:
                    if (H_i.induced_substructure(comb)).is_isomorphic(P_i):
                        p_sum += 1
                c_i += p_sum/len(H_combos)
            c_values.append(c_i)

        variables = ['p' + str(n) for n in range(1,num_graphs+1)]
        c_i_variables = [p if c_values[variables.index(p)] == 1 
                         else str(c_values[variables.index(p)]) + p 
                         if c_values[variables.index(p)] != 0 else 0 for p in variables]
        minimisation_list = [x for x in c_i_variables if x != 0]
        PSD_matrix = [[] for x in range(mat_dim)]

        for i in range(mat_dim):
            for j in range(mat_dim):
                entries_zeros = [variables[k] if H_matrices[k][i,j] == 1 
                                 else str(H_matrices[k][i,j]) + variables[k] if H_matrices[k][i,j] != 0 
                                 else 0 for k in range(num_graphs)]
                entries = [entry for entry in entries_zeros if entry != 0]
                PSD_matrix[i].extend([(' + '.join(map(str, entries))) if (' + '.join(map(str, entries))) != '' else 0])

        # Save written output to file:
        with open(savedirectory + "output.txt", "a") as f:
            print('\nIn order to calculate the inducibility, we solve the following semidefinite program:\n' + 
                  '\nVariables: ' + ', '.join(map(str, variables)) + '\nObjective: maximise ' + 
                  ' + '.join(map(str, minimisation_list)) + '\nConstraints:\n' + ', '.join(map(str, variables)) + 
                  ' >= 0' + '\n' + ' + '.join(map(str, variables)) + ' = 1' + '\n' + 
                  '[' + '\n'.join(map(str,PSD_matrix)) + '] >= 0\n',file = f)


        # SDP output:
        # This section uses strings and exec, which is not optimal but there was no functionality for this
        constraint_string = 'p.add_constraint('
        objective_string = ['c_values[' + str(i) + ']*x[' + str(i) + ']' for i in range(num_graphs)]
        sum_strings = ['matrix([[1]])*x[' + str(i) + ']' for i in range(num_graphs)]
        matrix_string = ['H_matrices[' + str(i) + ']*x[' + str(i) + ']' for i in range(num_graphs)]

        # Statement of SDP:
        p = SemidefiniteProgram(solver='cvxopt',maximization=True)
        x = p.new_variable()

        # Objective:
        exec('p.set_objective(' + ' + '.join(objective_string) + ')')

        # Constraints:
        exec(constraint_string + ' + '.join(sum_strings) + ' == matrix([[1]]))')

        for i in range(len(k_graphs_list)):
            p.add_constraint(matrix([[1]])*x[i] >= 0)

        exec(constraint_string + ' + '.join(matrix_string) + ' >= matrix.zero(len(all_flags)))')

        # Solve SDP and output solution to file:
        with open(savedirectory + "output.txt", "a") as f:
            print('This gives an optimal value of: ' + str(p.solve()) + '\n',file = f)
            print('The corresponding variable values are: ' + str(p.get_values(x)),file = f)


    # Turan density problem
    if problem == 'D':
        
        # Written output:
        w = sum([i for i in range(1,mat_dim+1)]) 
        H_coefficients = [[] for x in range(num_graphs)]
        q_variables = [None] * w
        Q = [[] for l in range(mat_dim)] 
        Q_matrices = []
        for s in range(mat_dim):
            for t in range(s,mat_dim):
                M = matrix.zero(mat_dim)
                M[s,t] = 1
                M[t,s] = 1
                Q_matrices.append(M)
                if s == t:
                    q_variables[mat_dim*s + t - sum(range(s+1))] = 'q' + str(s) + str(t)
                    Q[s].append('q' + str(s) + str(t))      
                else:
                    q_variables[mat_dim*s + t - sum(range(s+1))] = 'q' + str(s) + str(t)
                    Q[s].append('q' + str(s) + str(t))
                    Q[t].append('q' + str(s) + str(t))

        # Densities d(H) for each graph in k_graphs_list
        dens_mat = matrix(RR,len(k_graphs_list))
        for H in k_graphs_list:
            dens_mat[k_graphs_list.index(H),k_graphs_list.index(H)] = H.num_blocks()/math.comb(k,r)

        H_constraints = [[] for x in range(num_graphs)]
        for i in range(num_graphs):
            H_coefficients[i] = sum([[H_matrices[i][s,t] if s == t 
                                      else H_matrices[i][s,t] + H_matrices[i][t,s] 
                                      for t in range(s,mat_dim)] for s in range(mat_dim)],[])
            constraint = [q_variables[r] if H_coefficients[i][r] == 1 
                          else str(H_coefficients[i][r]) + q_variables[r] 
                          if H_coefficients[i][r] != 0 else 0 for r in range(w)]
            H_constraints[i] = [c for c in constraint if c != 0] 
            if dens_mat[i,i] != 0:
                H_constraints[i].append(str(dens_mat[i,i]))


        # Save written output to file:
        with open(savedirectory + "output.txt", "a") as f:
            print('\nIn order to calculate the Turan density, we solve the following semidefinite program:\n' + 
            '\nVariables: ' + ', '.join(q_variables) + ', k' + 
            '\nObjective: minimise k' + '\nConstraints: ',file = f)
            for i in range(num_graphs):
                print(' + '.join(H_constraints[i]) + ' <= k',file = f)
            print('[' + '\n'.join(map(str,Q)) + '] >= 0\n',file = f) 


        # SDP output:
        # This section uses strings and exec, which is not optimal but there was no functionality for this

        # Statement of SDP:
        p = SemidefiniteProgram(maximization=False)
        q = p.new_variable()

        # Objective:
        p.set_objective(q[w])

        # Constraints PSD matrix:
        constraint_matrices = []
        constraint_string = 'p.add_constraint('
        Q_matrix_string_list = []
        constraint_string_list = []
        for i in range(w):
            c_mat = matrix(RR,num_graphs,num_graphs)
            for j in range(num_graphs):
                c_mat[j,j] = H_coefficients[j][i]
            constraint_matrices.append(c_mat)
            Q_matrix_string_list.append('q[' + str(i) + ']*Q_matrices[' + str(i) + ']')
            constraint_string_list.append('q[' + str(i) + ']*constraint_matrices[' + str(i) + ']')

        # Constraints:
        exec(constraint_string + 'dens_mat + ' + ' + '.join(constraint_string_list) + ' <= matrix.identity(len(k_graphs_list))*q[w])')
        exec(constraint_string + ' + '.join(Q_matrix_string_list) + ' >= 0)')

        # Solve SDP and output solution to file:
        with open(savedirectory + "output.txt", "a") as f:
            print('This gives an optimal value of: ' + str(p.solve()) + '\n',file = f)
            print('The corresponding variable values are: ' + str(p.get_values(q)),file = f)

    return('')