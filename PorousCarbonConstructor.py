import streamlit as st
import random as ran
import numpy as np
#import scipy
#from scipy import stats as st
import time
#from periodictable import H,C 

st.title("Porous Carbon Constructor")

st.sidebar.markdown("## Paramater initialization")

## Initialize Data
st.sidebar.number_input("Number of Atoms:", min_value=500, step = 50, key="num_atoms")
st.sidebar.number_input("Foam Density:", min_value=0.05, value = 0.5, step = 0.05, key="density")
st.sidebar.number_input("Number of Pores:", min_value=1, value = 25, step=1, key="num_pores")
st.sidebar.slider("Porosity:", min_value=0.0,max_value=1.0, value =0.5, key='porosity')
st.sidebar.slider("Maximum Pore Size (e.g 0.5 will be half the box lenght):", min_value=0.01,max_value=0.99, value =0.5, key='max_pore_size_value')
st.sidebar.slider("Pore Overlap:", min_value=0.00,max_value=0.99, value =0.3, key='pore_overlap')
st.sidebar.slider("C-C Cutoff", min_value=1.0,max_value=1.4, step = 0.1, value =1.2, key='cutoff')

num_atoms_arr = np.array([st.session_state.num_atoms], dtype = np.int64)

### Naming convention for saved files ##########################################################
stringDensity = str(st.session_state.density).replace(".","p")+"gcc_"
stringNumAtoms = str(st.session_state.num_atoms)+"atoms_"
stringFoamOverlap = str(st.session_state.pore_overlap).replace(".","p")+"overlap"
stringNoFoam = str(st.session_state.num_pores)+"pores_"

atoms_vasp = "POSCAR_"+stringNumAtoms+stringDensity+stringNoFoam+stringFoamOverlap
pores_vasp = "POSCAR_PORES_"+stringNumAtoms+stringDensity+stringNoFoam+stringFoamOverlap
atoms_and_pores_xyz = "ovito_"+stringNumAtoms+stringDensity+stringNoFoam+stringFoamOverlap+".xyz"
#################################################################################################


# You can access the value at any point with:
st.session_state.num_atoms
st.session_state.density
st.session_state.num_pores
st.session_state.porosity
st.session_state.max_pore_size_value
st.session_state.pore_overlap
st.session_state.cutoff



##################### FUNCTIONS ##############################################################################################

def boxSize(*,density=st.session_state.density):
    '''
    This function predicts the box size for the porous carbon model, in units of angstrom
    '''

    atom_mass_amu = np.array([12.0107])
    atom_mass_gram = [atom*1.66054e-24 for atom in atom_mass_amu]  #Convert amu to grams (1 amu = 1.66054e-24 g)

    total_mass = sum([atom_mass_gram[i]*num_atoms_arr[i] for i in range(np.size(atom_mass_gram))])
    volume = total_mass / density

    box_cm = volume**(1/3) #(cm)
    ##convert box lenght  in cm to armstrong
    box_arm = box_cm/1e-8
    st.write(f"The box lenght is {box_arm:.2f} \u212B")
    return box_arm, volume


def uniformPore(*, max_pore_size=st.session_state.max_pore_size_value):
    '''
    This function obtains a uniform pore distribution
    The max_pore_size paramter can be changed to desired maximum pore size
    '''
    return np.random.rand() * (box*max_pore_size)


def chiPore(*, max_pore_size = st.session_state.max_pore_size_value , df = 1,loc = 0, scale = .25, size=1):
    '''
    This function obtains a chi distributed pore distribution (https://en.wikipedia.org/wiki/Chi_distribution)
    The default parameter aguments are set to resemble a half-Gaussian
    The scale of .25 is set so that the maximum is around 1, to get the max allowed pore sized
    '''
    return np.random.chisquare(df, size)[0] * (box*max_pore_size)


def betaPore(*, max_pore_size=st.session_state.max_pore_size_value , a = 4, b = 4):
    '''
    This function obatins a beta distribution (https://en.wikipedia.org/wiki/Beta_distribution)
    The parameters, a and b, set at 4 give a non-negative bell-shaped distribution, akin to a Gaussian Distribution
    '''
    return np.random.beta(4,4) * (box*max_pore_size)

def poreCreator(*, porosity = st.session_state.porosity, num_pores=st.session_state.num_pores, epsilon=0.1, pore_dist_kind = 3):
    '''
    This function creates the pores in a desired distribution and porosity
    
    porosity is the desired porosity, example; 0.5
    num_pores is the number of pores you want to sample
    epsilon ensures that we select a distribution that is less the desired porosity. We set this to 0.1
    pore_dist_kind specify what pore distribution is preffered example 3 (default) is the beta distribution
    '''

    poreRadii_list = []
    poreVolume_sum = 0
    i = 0           #iterator for our loop
    restarts = 0    #tells how many times this function is restarted when the distribution surpases the desired porsity
    porosity_threshold = porosity - epsilon


    while i < num_pores:
    
        if poreVolume_sum <= porosity_threshold:
            
            if pore_dist_kind == 1:
                name = "Uniform Distribution"
                pore_radius = uniformPores()
                
            elif pore_dist_kind == 2:
                name = "Chi Distribution"
                pore_radius = chiPore()
                
            else:
                name = "Beta Distribution"
                pore_radius = betaPore()
            
            poreVolume_sum += pore_radius**3/box**3
            poreRadii_list.append(pore_radius)
            i += 1
        
              
            if poreVolume_sum > porosity:
                restarts += 1
                print(f"Number of restarts = {restarts}", end = '\r')
                poreVolume_sum = 0
                i = 0
                poreRadii_list = []
            
            
            # Ensure the last pore sums to the desired porosity
            if i == num_pores-1:
                print(f"Sampled up to {i} pores. Current Porosity is: {poreVolume_sum:.2f}")
                last_poreVolume = abs((porosity - poreVolume_sum)) * box**3
                print(f"Added porosity fraction is {(last_poreVolume/box**3):.2f}")
                poreRadii_list.append(last_poreVolume**(1/3))
                poreVolume_sum += last_poreVolume/box**3
                print(f"Added last pore. Final porosity is: {poreVolume_sum:.2f}")
                break
            
   
        if poreVolume_sum > porosity_threshold:
            restarts += 1
            print(f"Number of restarts = {restarts}", end='\r')
            poreVolume_sum = 0
            i = 0
            poreRadii_list = []
        
    print()
    print(f"                 *** Pore Radii List for {len(poreRadii_list)} Pores ***")
    print(np.round(poreRadii_list,2))
    
    return name, poreRadii_list
##########################################################################################################

# Input parameters for boxSize
st.header('Box Size Calculation')

if st.button('Calculate Box Size'):
    box, volume = boxSize(density=st.session_state.density)
    st.write(f"The box length is {box:.2f} Å")
    st.write(f"The volume is {volume:.2e} cm³")

# Input parameters for pore distribution
st.header('Pore Distribution')
pore_dist_kind = st.selectbox('Pore Distribution Kind', options=[1, 2, 3], format_func=lambda x: ["Uniform", "Chi", "Beta"][x-1])

if st.button('Create Pores'):
    name, poreRadii_list = poreCreator(pore_dist_kind=pore_dist_kind)
    st.write(f"Pore Distribution: {name}")
    st.write(f"Pore Radii List: {np.round(poreRadii_list, 2)}")
