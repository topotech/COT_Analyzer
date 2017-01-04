# *-* coding: utf-8 *-*
"""
COTlib
======

Python package for complex systems analysis through Chemical Organization
Theory.

This package provides three classes needed to build a reaction network:
    1.- species:    Species
    2.- reaction:   Reactions
    3.- renet:      Reaction Networks

This package also provides two global sets containing all the species
and reactions created, called M and R respectively.

"""

import cvxopt as co
import numpy as np
import re
import hasselib
import networkx as nx
import matplotlib.pyplot as plt
import os

"""
Global variables and functions
==============================
"""

# Global sets
M = list()
R = list()

# Display all elements in M
def show_M():
    for i, elem in enumerate(M):
        print "["+str(i)+"]: "+str(elem)

# Display all reactions in R
def show_R():
    for i,reac in enumerate(R):
        print "["+str(i)+"]: "+str(reac)

# Remove species from M
def rm_species(name):
    if not(isinstance(name,str)):
        raise TypeError("rm_species() expects a string")
    for elem in M:
        if name == elem.name:
            M.remove(elem)
            return True
    return False

# Verify of species exists in M: return boolean and index (-1 if false)
def in_M(name):
    if not(isinstance(name,str)):
        raise TypeError("in_M() expects a string")
    for i,elem in enumerate(M):
        if name == elem.name:
            return True, i
    return False, -1

def get_stats_of_lens():
    lens = np.array([])
    for r in R:
        lens=np.append(len(set(r.reac_list).union(set(r.prod_list))),lens)
    return np.std(lens), np.mean(lens), np.median(lens)

#Remove all reactions where number of species participating is greater than p*std+median (or mean)
def filter_by_std(p,use_mean=False):
    std,mean,median = get_stats_of_lens()
    new_R = []
    to_remove = []
    for r in R:
        if use_mean:
            #if len(set(r.reac_list).union(set(r.prod_list))) <= p*std+mean:
            if len(set(r.reac_list)) <= p*std+mean:
                new_R.append(r)
            else:
                to_remove.append(r)
        else:
            #if len(set(r.reac_list).union(set(r.prod_list))) <= p*std+median:
            if len(set(r.reac_list)) <= p*std+median:
                new_R.append(r)
            else:
                to_remove.append(r)
    del R[:]
    R.extend(new_R)
    for s in M:
        for r in to_remove:
            try:
                s.myReactions.remove(r)
            except ValueError:
                pass


#Remove all reactions where number of species participating is greater than tol+median (or mean)
def filter_by_absolute(tol,use_mean=False):
    std,mean,median = get_stats_of_lens()
    new_R = []
    to_remove = []
    for r in R:
        if use_mean:
            #if len(set(r.reac_list).union(set(r.prod_list))) <= tol+mean:
            if len(set(r.reac_list)) <= tol+mean:
                new_R.append(r)
            else:
                to_remove.append(r)
        else:
            #if len(set(r.reac_list).union(set(r.prod_list))) <= tol+median:
            if len(set(r.reac_list)) <= tol+median:
                new_R.append(r)
            else:
                to_remove.append(r)
    del R[:]
    R.extend(new_R)
    for s in M:
        for r in to_remove:
            try:
                s.myReactions.remove(r)
            except ValueError:
                pass

def get_inflow_indexes():
    result = []
    for i,r in enumerate(R):
        if len(r.reactants) == 0:
            result.append(i)
    return result

"""
Classes and their methods
=========================
"""


# Class species
class species:
    def __init__(self, name):
        if not(isinstance(name,str)):
            raise TypeError("species constructor expects a string")

        for elem in M:
            if name == elem.name:
                raise ValueError("there can't be two elements with the same name")
        self.myReactions = list()
        self.name = name
        M.append(self)

    def __repr__(self):
        return self.name

    def ins_reaction(self,myreaction):
        if not(isinstance(myreaction,reaction)):
            raise TypeError("ins_reaction() expects an instance of reaction class")
        self.myReactions.append(myreaction)


# Class reaction
class reaction:

    def __init__(self, input):
        if not(isinstance(input,str) or isinstance(input,reaction_input)):
            raise TypeError("reaction constructor expects a string or an instance of \'reaction_input\' ")

        if isinstance(input,str):

            self.repr = input.replace(" ", "")

            pattern = re.compile(r"(.*)(?:->|→|⇒)(.*)")
            match = re.match(pattern, self.repr)

            if not(match):
                raise ValueError("The reaction string was badly formatted")
            if len(match.groups()) != 2:
                raise ValueError("The reaction string was badly formatted")

            str_reactants = match.groups()[0].split("+")
            str_products  = match.groups()[1].split("+")

            self.reactants = set()
            self.products   = set()
            pattern = re.compile(r"([0-9]*)(\w*)")
            #pattern = re.compile(r"(\d*\.?\d*)(\w*)")


            if len(match.groups()[0]) > 0:
                for string in str_reactants:
                    match2 = re.match(pattern, string)
                    if not(match2):
                        raise ValueError("The reactants were badly formatted")
                    if len(match2.groups()[1]) <= 0:
                        raise ValueError("The reactants were badly formatted")

                    if len(match2.groups()[0]) == 0:
                        self.reactants.add((1,match2.groups()[1]))
                    else:
                        self.reactants.add((int(match2.groups()[0]),match2.groups()[1]))
                        #self.reactants.add((float(match.groups()[0]),match.groups()[1]))

            if len(match.groups()[1]) > 0:
                for string in str_products:
                    match = re.match(pattern, string)
                    if not(match):
                        raise ValueError("The products were badly formatted")
                    if len(match.groups()[1]) <= 0:
                        raise ValueError("The products were badly formatted")

                    if len(match.groups()[0]) == 0:
                        self.products.add((1,match.groups()[1]))
                    else:
                        self.products.add((int(match.groups()[0]),match.groups()[1]))
                        #self.products.add((float(match.groups()[0]),match.groups()[1]))

            self.reac_list = [item[1] for item in self.reactants]
            self.prod_list = [item[1] for item in self.products]
            self.reac_list.sort(key=lambda x:in_M(x)[1])
            self.prod_list.sort(key=lambda x:in_M(x)[1])
            if len(self.reac_list) != len(set(self.reac_list)):
                raise ValueError("The reactants were badly formatted: there are repeated species")

            if len(self.prod_list) != len(set(self.prod_list)):
                raise ValueError("The products were badly formatted: there are repeated species")

            for reac in R:
                if self.reactants == reac.reactants and self.products == reac.products:
                    raise ValueError("there can't be two equal reactions")

            for tup in self.reactants.union(self.products):
                boolean,index = in_M(tup[1])
                if boolean:
                    M[index].ins_reaction(self)
                else:
                    species(tup[1])
                    boolean,index = in_M(tup[1])
                    M[index].ins_reaction(self)

            R.append(self)
        elif isinstance(input,reaction_input):
            self.repr = input.repr


            self.reactants = input.reactants
            self.products= input.products
            self.reac_list = [item[1] for item in self.reactants]
            self.reac_list = [item[1] for item in self.reactants]
            self.prod_list = [item[1] for item in self.products]
            self.reac_list.sort(key=lambda x:in_M(x)[1])
            self.prod_list.sort(key=lambda x:in_M(x)[1])

            if len(self.reac_list) != len(set(self.reac_list)):
                raise ValueError("The reactants were badly formatted: there are repeated species")

            if len(self.prod_list) != len(set(self.prod_list)):
                raise ValueError("The products were badly formatted: there are repeated species")

            for reac in R:
                if self.reactants == reac.reactants and self.products == reac.products:
                    raise ValueError("there can't be two equal reactions")

            for tup in self.reactants.union(self.products):
                boolean,index = in_M(tup[1])
                if boolean:
                    M[index].ins_reaction(self)
                else:
                    species(tup[1])
                    boolean,index = in_M(tup[1])
                    M[index].ins_reaction(self)

            R.append(self)


    def __repr__(self):
        return self.repr

    def __eq__(self,other):
        return self.reactants == other.reactants and self.products == other.products
    def __hash__(self):
        return hash((tuple(self.reactants),tuple(self.products)))


# Class reaction network
class renet:
    S_x = np.array([])

    def __init__(self, indexes):
        if not(isinstance(indexes,list)):
            raise TypeError("renet constructor expects a list of integers")
        if not (all(isinstance(x, int) for x in indexes)):
            raise TypeError("renet constructor expects a list of integers")

        self.X = list()
        self.X_names = set()
        self.R_x = list()

        for i in indexes:
            self.X.append(M[i])
            self.X_names.add(M[i].name)

        '''
        for elem in self.X:
            for reaction in elem.myReactions:
                self.ins_reaction(reaction)

            self.cln_Rx()
        '''

        for elem in self.X:
            for reaction in elem.myReactions:
                    self.ins_reaction(reaction)

        self.createA()
        self.createB()

        self.S = self.B - self.A


    # This let the user to add reactions with empty reactants outside X
    def add_species_for_closing(self, indexes):
        if not(isinstance(indexes,list)):
            raise TypeError("add_species_for_closing() expects a list of integers")
        if not (all(isinstance(x, int) for x in indexes)):
            raise TypeError("add_species_for_closing() expects a list of integers")

        for i in indexes:
            if not (M[i].name in self.X_names):
                self.X.append(M[i])
                self.X_names.add(M[i].name)
        '''
        for elem in self.X:
            for reaction in elem.myReactions:
                self.ins_reaction(reaction)

        self.cln_Rx()
        '''
        for elem in self.X:
            for reaction in elem.myReactions:
                    self.ins_reaction(reaction)


        #self.createA()
        #self.createB()

        #self.S = self.B - self.A


    # This let the user to add reactions with empty reactants outside X
    def add_inflow(self, indexes):
        if not(isinstance(indexes,list)):
            raise TypeError("add_inflow() expects a list of integers")
        if not (all(isinstance(x, int) for x in indexes)):
            raise TypeError("add_inflow() expects a list of integers")


        for i in indexes:
            if len(R[i].reactants) == 0:
                self.ins_reaction(R[i])
            else:
                raise ValueError("Malformed input. All reactant sets need to be empty.")

        #self.cln_Rx()
        self.createA()
        self.createB()

        self.S = self.B - self.A

    #cln_Rx(): Clean reaction: just leave the reactions which can be fired by X
    def cln_Rx(self):
        aux = list()
        for myreaction in self.R_x:
            reactants = set([item[1] for item in myreaction.reactants])
            if not (reactants.issubset(self.X_names)):
                aux.append(myreaction)
        for item in aux:
            self.R_x.remove(item)


    def ins_reaction(self,r):
        if not set(self.X_names).issuperset(set(r.reac_list)):
            return False
        for reac in self.R_x:
            if (r.reactants == reac.reactants) and (r.products == reac.products):
                return False
        self.R_x.append(r)
        return True

    #createA(): Create reactants matrix
    def createA(self):

        self.A = np.zeros(  (len(self.X) , len(self.R_x) ) )

        for i,r in enumerate(self.R_x):

            for j,e in enumerate(self.X):

                for m in r.reactants:
                    if m[1] == e.name:
                        self.A[j][i] = m[0]


    #createB(): Create products matrix
    def createB(self):

        self.B = np.zeros(  (len(self.X) , len(self.R_x) ) )

        for i,r in enumerate(self.R_x):

            for j,e in enumerate(self.X):

                for m in r.products:
                    if m[1] == e.name:
                        self.B[j][i] = m[0]



    #close(): Turn X into it's closure (add all the species produced by R_x into X).
    def close(self):
        #print "Hola"
        counter = 0
        checked_reactions = set()
        while not(self.isClosed()):
            counter += 1
            indexes = list()
            producedSet = set()

            for r in self.R_x:
                if not (r in checked_reactions):
                    checked_reactions.add(r)
                    aux = set([item[1] for item in r.products])
                    producedSet = producedSet.union(aux)

            for name in producedSet:
                indexes.append( in_M(name)[1] )
            '''
            #Reinitializing the reaction network
            for i in indexes:
                if not (M[i].name in self.X_names):
                    self.X.append(M[i])
                    self.X_names.add(M[i].name)
            '''
            self.add_species_for_closing(indexes)


            '''
            for elem in self.X:
                for reaction in elem.myReactions:
                    self.ins_reaction(reaction)
            '''

            #self.cln_Rx() #CUELLO DE BOTELLA
        #
        # print "chao"
        self.createA()
        self.createB()

        self.S = self.B - self.A
        return counter



    #show_x(): Display a list with all elements in X and its indexes
    def show_X(self):
        for i, elem in enumerate(self.X):
            print "["+str(i)+"]: "+str(elem)

    #show_Rx(): Display all reactions in R_x
    def show_Rx(self):
        for i,reac in enumerate(self.R_x):
            print "["+str(i)+"]: "+str(reac)

    #isSM(): Verify if this reaction network is self-mantaining
    def isSM(self):
        #if X is empty, then it's SM
        if len(self.X) == 0:
            return True

        #vars is the numbers of columns of S. Represents
        # the  number of components of the unknow flux vector 'v'
        vars = np.shape(self.S)[1]

        #eqs1 is the numbers of rows of S. Represents
        # the number of inequations introduced in the LP problem.
        eqs1 = np.shape(self.S)[0]

        #eqs2 is the numbers of columns of S. Represents
        # the lower bound inequations introduced in the LP problem.
        eqs2 = vars

        I = np.identity(vars)
        A = self.S
        A = np.concatenate((A, I), axis=0)
        A = -A
        b0 =  np.zeros(eqs1)
        b1 = -np.ones(eqs2)
        b = np.concatenate((b0, b1), axis=0)

        c = np.ones(vars)

        A = co.matrix(A)
        b = co.matrix(b)
        c = co.matrix(c)

        co.solvers.options['show_progress']=False
        sol=co.solvers.lp(c,A,b)

        return bool(sol['x'])



    #isOverproduced(i): Verify if species m_i is overproduced. Use show_X() to list all indexes.
    def isOverproduced(self,i):

        if i >= len(self.X) or i < 0:
            raise ValueError("Index out of bound. A non-negative integer less than len(X) expected.")

        #vars is the numbers of columns of S. Represents
        # the  number of components of the unknow flux vector 'v'
        vars = np.shape(self.S)[1]

        #eqs1 is the numbers of rows of S. Represents
        # the number of inequations introduced in the LP problem.
        eqs1 = np.shape(self.S)[0]

        #eqs2 is the numbers of columns of S. Represents
        # the lower bound inequations introduced in the LP problem.
        eqs2 = vars

        I = np.identity(vars)
        A = self.S
        A = np.concatenate((A, I), axis=0)
        A = -A

        b =  np.zeros(eqs1+eqs2)
        b[i] = -1

        c = np.ones(vars)

        A = co.matrix(A)
        b = co.matrix(b)
        c = co.matrix(c)

        co.solvers.options['show_progress']=False
        sol=co.solvers.lp(c,A,b)

        return bool(sol['x'])

    #show_x(): Display a list with all overproduced elements in X
    def showOverproduced(self):
        boolean = True
        for i, elem in enumerate(self.X):
            if(self.isOverproduced(i)):
                print "["+str(i)+"]: "+str(elem)
                boolean = False

        if (boolean): print "None"


    #isClosed(): Return True if X is closed (False if not).
    def isClosed(self):
        producedSet = set()
        for r in self.R_x:
            aux = set([item[1] for item in r.products])
            producedSet = producedSet.union(aux)
        #if the elements produced by R_x is subset of X
        return producedSet.issubset(set(self.X_names))

    #isSSM(): Return True if X is semi self-maintaining (False if not).
    def isSSM(self):
        producedSet = set()
        consumedSet = set()
        for r in self.R_x:
            aux1 = set([item[1] for item in r.products])
            producedSet = producedSet.union(aux1)
            aux2 = set([item[1] for item in r.reactants])
            consumedSet = consumedSet.union(aux2)


        #if X is subset of the elements produced by R_x
        return consumedSet.issubset(producedSet)

    #Obtain a hierchy of closed sets (returns a Networkx graph)
    def reactant_closure (self, add_separables=False, verbose=False, render=False,save_plot=False,save_graph=False,print_log=False):
        total_of_species = len(self.X)
        # Save global indexes to inflow (empty) reactions.
        R_inflow = list()

        loiReactants = set() #set of tuples containing the indexes of each reaction

        loiReactants.add( () ) #Including the empty set
        for r in self.R_x:
            loiReactants.add(tuple([in_M(item)[1] for item in r.reac_list]))


        #Intializing "tree":
        for i,my_reac in enumerate(self.R_x):
            if len(my_reac.reactants) == 0:
                for j,global_reac in enumerate(R):
                    if my_reac.reactants == global_reac.reactants and my_reac.products == global_reac.products:
                        R_inflow.append(j)
        tree_of_closed_sets = list()
        if print_log:
            i = 0
            while os.path.exists("graph_log"+str(i)+".log"):
                i += 1
            my_log_filename = "graph_log"+str(i)+".log"
            my_log_string = ""


        for S in loiReactants:
            l = list(S)
            flag = False
            curr_set = renet(l)
            aux= set([in_M(item)[1] for item in curr_set.X_names]) #Add initials
            tree_of_closed_sets.append(node_closure_tree(aux,0))
            if print_log:
                my_log_string += self.stringify_node(tree_of_closed_sets[-1])


            curr_index = len(tree_of_closed_sets)-1

            curr_set.add_inflow(R_inflow)
            curr_set.close()
            aux_closed= set([in_M(item)[1] for item in curr_set.X_names]) #Add closed

            for element in reversed(tree_of_closed_sets):
                if element.my_set == aux_closed:
                    flag = True
            if flag == False:
                tree_of_closed_sets[curr_index].closure.append(aux_closed)
                tree_of_closed_sets.append(node_closure_tree(aux_closed,1))
                if print_log:
                    my_log_string += self.stringify_node(tree_of_closed_sets[-1])
                tree_of_closed_sets[curr_index+1].closure.append(aux_closed)
            else:
                del curr_set
                tree_of_closed_sets[curr_index].closure.append(aux_closed)
                if aux_closed == aux:
                    tree_of_closed_sets[curr_index].level = 1
                continue
            del curr_set

        tree_of_closed_sets.sort(key=lambda x: x.level)
        #Calculate the range inside the list containing the level 0 elements.
        c=0
        while (tree_of_closed_sets[c].level < 1):
            c+=1
        root_init = c
        while (tree_of_closed_sets[c].level == 1):
            c+=1
            if c == len(tree_of_closed_sets): break
        root_final = c
        closed_reactants = tree_of_closed_sets[root_init:root_final]

        c=0
        k=1
        if verbose:
            max_number_of_species = max([len(node.my_set) for node in closed_reactants])

        while (k <= tree_of_closed_sets[-1].level):
            if verbose:
                print "Starting with level " + str(k)
            if print_log:
                my_log_string += "Starting with level " + str(k)+"\n"

            while (tree_of_closed_sets[c].level < k):
                c+=1
            init = c
            while (tree_of_closed_sets[c].level == k):
                c+=1
                if c == len(tree_of_closed_sets): break
            final = c

            l = tree_of_closed_sets[init:final]
            if verbose:
                self.verbose_reactantclosure((float(max_number_of_species)/total_of_species),len(tree_of_closed_sets))
            if print_log:
                self.log_printer(my_log_filename,my_log_string)
                my_log_string = ""


            for i in l:
                if verbose:
                    self.verbose_reactantclosure(float(max_number_of_species)/total_of_species ,len(tree_of_closed_sets))
                if print_log:
                    self.log_printer(my_log_filename,my_log_string)
                    my_log_string = ""
                for j in closed_reactants:
                    if i.my_set.issuperset(j.my_set): continue

                    union_set = i.my_set.union(j.my_set)

                    if add_separables == False:
                        U_renet = renet(list(union_set))
                        A_renet = renet(list(j.my_set))
                        B_renet = renet(list(i.my_set))
                        U_renet.add_inflow(R_inflow)
                        A_renet.add_inflow(R_inflow)
                        B_renet.add_inflow(R_inflow)
                        U_reactions = set([str(reac) for reac in U_renet.R_x])
                        A_reactions = set([str(reac) for reac in A_renet.R_x])
                        B_reactions = set([str(reac) for reac in B_renet.R_x])
                        del U_renet; del A_renet; del B_renet
                        Uisseparable = ( U_reactions == A_reactions.union(B_reactions) )

                    else: Uisseparable = False

                    if not Uisseparable:
                        flag = False
                        '''
                        for element in tree_of_closed_sets:
                            if element.my_set == union_set:
                                flag = True
                                break
                        '''
                        if flag == False:
                            '''
                            flag = False
                            for element in i.united_to_form:
                                if element == union_set:
                                    flag = True
                                    break
                            if flag == False: i.united_to_form.append(union_set)

                            flag = False
                            for element in j.united_to_form:
                                if element == union_set:
                                    flag = True
                                    break
                            if flag == False: j.united_to_form.append(union_set)

                            tree_of_closed_sets.append(node_closure_tree(union_set,k+1)) #Add union to queue of next iteration
                            if print_log:
                                my_log_string += self.stringify_node(tree_of_closed_sets[-1])
                            '''
                            curr_set = renet(list(union_set))
                            curr_set.add_inflow(R_inflow)
                            curr_set.close()
                            closed_union_set = set([in_M(item)[1] for item in curr_set.X_names]) #Get set of indexes from closed
                            del curr_set
                            #if closed_union_set != union_set:
                            if True:
                                #tree_of_closed_sets[-1].closure.append(closed_union_set)
                                flag = False
                                for element in reversed(tree_of_closed_sets):
                                    if element.my_set == closed_union_set:
                                        flag = True
                                        break
                                if flag == False:
                                    #tree_of_closed_sets.pop() #Remove last added node and...
                                    tree_of_closed_sets.append(node_closure_tree(closed_union_set,k+1)) #...replace it by its closure
                                    if verbose:
                                        candidate_max_number_of_species = len(closed_union_set)
                                        if max_number_of_species < candidate_max_number_of_species:
                                            max_number_of_species = candidate_max_number_of_species

                                    tree_of_closed_sets[-1].closure.append(closed_union_set)
                                    if print_log:
                                        my_log_string += self.stringify_node(tree_of_closed_sets[-1])
                            #else:
                            #    tree_of_closed_sets[-1].closure.append(closed_union_set)
                        '''
                        else:
                            flag = False
                            for element in i.united_to_form:
                                if element == union_set:
                                    flag = True
                                    break
                            if flag == False: i.united_to_form.append(union_set)

                            flag = False
                            for element in j.united_to_form:
                                if element == union_set:
                                    flag = True
                                    break
                            if flag == False: j.united_to_form.append(union_set)
                            continue
                        '''
            tree_of_closed_sets.sort(key=lambda x: x.level)

            if verbose:
                self.verbose_reactantclosure(float(max_number_of_species)/total_of_species ,len(tree_of_closed_sets))
            if print_log:
                self.log_printer(my_log_filename,my_log_string)
                my_log_string = ""
            if verbose:
                print "Ending with level " + str(k)
            k += 1


        #Note to self: Verify which sets are closed during the loop (we're overkilling
        # it if we have to iterate the whole set again).
        for element in tree_of_closed_sets:
            if len(element.closure) == 0:
                element.IsClosed = True
                print "OJO"
                continue
            if element.closure[0] == element.my_set:
                element.IsClosed = True
            else:
                element.IsClosed = False


        '''
        for i in tree_of_closed_sets:
            print i.my_set
            if len(i.united_to_form)>0:
                print "  Uniones:"
                for j in i.united_to_form:
                    print "     "+str(j)
            else:
                print "  Uniones:"
                print "     set([])"
            if len(i.closure)>0:
                print "  closure:"
                for j in i.closure:
                    print "     "+str(j)
            else:
                print "  closure:"
                print "     set([])"
            print "  closure:"
            print "     "+str(i.IsClosed)
        '''


        #Graph construction: Construction of the graph following algorithm beheaviour
        '''
        #G = nx.Graph()
        G = nx.DiGraph()
        sets_tree_of_closed_sets = [i.my_set for i in tree_of_closed_sets]
        levels = list()
        pos = dict()
        level_counter = list()
        closed_only = list()

        #Adding edges
        for i in tree_of_closed_sets:
            if i.IsClosed == True:
                levels.append(i.level)
                closed_only.append(i)
                if len(i.united_to_form) != 0:
                    for j in i.united_to_form:
                        index = sets_tree_of_closed_sets.index(j)
                        if tree_of_closed_sets[index].IsClosed:
                            G.add_edge(str(list(i.my_set)),str(list(j)))
                        else:
                            G.add_edge(str(list(i.my_set)),str(list(tree_of_closed_sets[index].closure[0])))
                else:
                    G.add_node(str(list(i.my_set)))

        #Building layout dictionary
        for i in range(k):
            level_counter.append(levels.count(i))
            if i == 0:
                init = 0
            else:
                init = final
            final = level_counter[-1] + init
            for j,element in enumerate(closed_only[init:final]):
                if level_counter[i] > 1:
                    if j%2 == 0:
                        pos[str(list(element.my_set))] = ( float(j)/(level_counter[i]-1) , float(i)/k )
                    else:
                        #pos[str(list(element.my_set))] = ( float(j)/(level_counter[i]-1) , float(i)/k )
                        pos[str(list(element.my_set))] = ( float(j)/(level_counter[i]-1) , float(i)/k + .08/k)
                else:
                    pos[str(list(element.my_set))] = ( 0.5 , float(i)/k )


        #nx.draw(G, pos, node_color='w', edge_color='r' , with_labels = True)
        #plt.show() # display
        '''


        #Adding edges
        closed_only = list()

        #Adding edges
        for i in tree_of_closed_sets:
            if i.IsClosed == True:
                closed_only.append(i)

        H = nx.Graph()
        closed_only.sort(key=lambda x: len(x.my_set))
        cardinality_counter =  [len(element.my_set) for element in closed_only]
        pos = dict()

        i = 0
        slice_points = [0]
        saved = cardinality_counter[0]
        while i < len(cardinality_counter):
             if cardinality_counter[i] != saved:
                     saved = cardinality_counter[i]
                     slice_points.append(i)
             i += 1
        slice_points.append(len(cardinality_counter))
        not_connected_yet = list()

        if  len(slice_points) != 2:
            for i in range(2,len(slice_points)):
                M = closed_only[slice_points[i-2]:slice_points[i-1]]
                N = closed_only[slice_points[i-1]:slice_points[i]]
                for sub in M:
                    flag = True #Test if conected
                    for sup in N:
                        if sub.my_set.issubset(sup.my_set):
                            H.add_edge(str(list(sub.my_set)),str(list(sup.my_set)))
                            flag = False
                    if flag != False:
                        not_connected_yet.append(sub)
                to_delete = list()
                for k,sub in enumerate(not_connected_yet):
                    flag = False
                    for sup in N:
                        if sub.my_set.issubset(sup.my_set):
                            H.add_edge(str(list(sub.my_set)),str(list(sup.my_set)))
                            flag = True
                    if flag == True:
                        to_delete.append(sub)

                while 0 < len(to_delete):
                    sub = to_delete[0]
                    not_connected_yet.remove(sub)
                    to_delete.remove(sub)

        else:
            for node in closed_only:
                H.add_node(str(list(node.my_set)))

        #Building layout

        #cardinalities = list(set(cardinality_counter))

        for i in range(len(slice_points)-1):
            M = closed_only[slice_points[i]:slice_points[i+1]]
            curr_cardinality = cardinality_counter[slice_points[i]]
            y = float(i)/(len(cardinality_counter)-1)
            for j in range(len(M)):
                if len(M) == 1:
                    x = 0.5
                else:
                    x = float(j)/(len(M)-1)
                pos[str(list(M[j].my_set))] = (x,y)

        if verbose:
            print "Closed sets found: " + str(H.number_of_nodes())
        if (render == True):
            nx.draw(H, pos, node_color='w', edge_color='r' , with_labels=True)
            plt.show() # display
        if (save_plot == True):
            nx.draw(H, pos, node_color='w', edge_color='r' , with_labels=True, node_size=60,font_size=8, figsize=(12,12))
            plt.savefig("graph.png")
        if (save_graph == True):
            nx.write_gml(H,"my.gml")

        return (H,tree_of_closed_sets)



    def verbose_reactantclosure(self,ratio, number_of_nodes):
        print "Curent ratio max_X_size/size_of_M  : "+( "%0.2f" % (ratio*100))
        print "Number of sets (nodes) in the graph: "+str(number_of_nodes)

    def stringify_node(self,node):
        string = node.string_for_log()

        string += "\n\n"
        return string

    def log_printer(self,log_filename,log_string):
        with open(log_filename, "a") as f:
             f.write(log_string)





class node_closure_tree:
    def __init__ (self,my_set,level):
        self.united_to_form = list()
        self.closure = list()
        self.my_set = my_set
        self.level = level
        self.IsClosed = None

    def __str__(self):
        string = str(list(self.my_set))+"\n"
        string += "united_to_form: "+str(self.united_to_form)+"\n"
        string += "Closure: "+str(self.closure)+"\n"
        string += "Level: "+str(self.level)
        return string

    def string_for_log(self):
        string = str(list(self.my_set))+"\n"
        string += "Level: "+str(self.level)
        return string

class reaction_input:
    def __init__(self):
        self.reactants = set()
        self.products = set()
        self.repr = ""

    def __repr__(self):
        return self.repr