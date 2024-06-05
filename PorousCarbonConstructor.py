import streamlit as st
import random as ran
import numpy as np
import pandas as pd
import time
import plotly
import plotly.express as px

st.title("Porous Carbon Constructor")
logo_url = "https://chinonsougwumadu.com/wp-content/uploads/2024/05/microsoftteams-image-17.jpg"
st.sidebar.image(logo_url)

# col3, col4 = st.columns(2)
# with col3:
#     st.sidebar.link_button("About Us", "https://daviddrabold.com/")
# with col4:
#     st.sidebar.link_button("Help", "https://daviddrabold.com/publications/")
st.sidebar.link_button("About Us", "https://daviddrabold.com/")
st.sidebar.link_button("Help", "https://daviddrabold.com/publications/")
st.sidebar.markdown("## Paramater initialization")




### To disable the generate key until the pore distribution button is clicked on
def disable(b):
    st.session_state["disabled"] = b


## Initialize Data
st.sidebar.number_input("Number of Atoms:", min_value=500, step = 50, key="num_atoms",on_change=disable, args=(True,), help="The number of C atoms required")
st.sidebar.number_input("Foam Density [g/cm$^3$]:", min_value=0.05, value = 0.5, step = 0.05, key="density",on_change=disable, args=(True,), help ="The desired foam density")
st.sidebar.number_input("Number of Pores:", min_value=1, value = 25, step=1, key="num_pores",on_change=disable, args=(True,), help="The number of pores required")
st.sidebar.slider("Porosity:", min_value=0.0,max_value=1.0, value =0.5, key='porosity',on_change=disable, args=(True,), help="Desired porosity of the model")
st.sidebar.slider("Maximum Pore Size", min_value=0.01,max_value=0.99, value =0.5, key='max_pore_size',on_change=disable, args=(True,), help="For example, 0.5 will be half the box lenght")
st.sidebar.slider("Pore Overlap:", min_value=0.00,max_value=0.99, value =0.3, key='pore_overlap',on_change=disable, args=(True,), help="Specicy if pore overlap is allowed")
st.sidebar.slider("Carbon Bonds Initial Cutoff [\u212B]", min_value=1.0,max_value=1.4, step = 0.1, value =1.2, key='cutoff',on_change=disable, args=(True,), help="C-C cutoff, 1.2 is ideal for optimal performance of app")


# You can access the value at any point with:
# st.session_state.num_atoms
# st.session_state.density
# st.session_state.num_pores
# st.session_state.porosity
# st.session_state.max_pore_size
# st.session_state.pore_overlap
# st.session_state.cutoff

st.session_state.num_atoms_arr = np.array([st.session_state.num_atoms], dtype = np.int64)

### Naming convention for saved files ##########################################################
stringDensity = str(np.round(st.session_state.density,2)).replace(".","p")+"gcc_"
stringNumAtoms = str(st.session_state.num_atoms)+"atoms_"
stringFoamOverlap = str(np.round(st.session_state.pore_overlap,2)).replace(".","p")+"overlap"
stringNoFoam = str(st.session_state.num_pores)+"pores_"

st.session_state.atoms_vasp = "POSCAR_"+stringNumAtoms+stringDensity+stringNoFoam+stringFoamOverlap
st.session_state.pores_vasp = "POSCAR_PORES_"+stringNumAtoms+stringDensity+stringNoFoam+stringFoamOverlap
st.session_state.atoms_and_pores_xyz = "ovito_"+stringNumAtoms+stringDensity+stringNoFoam+stringFoamOverlap+".xyz"
#################################################################################################

##################### FUNCTIONS ##############################################################################################

#def boxSize(*,density=st.session_state.density):
def boxSize(density):
    '''
    This function predicts the box size for the porous carbon model, in units of angstrom
    '''

    atom_mass_amu = np.array([12.0107])
    atom_mass_gram = [atom*1.66054e-24 for atom in atom_mass_amu]  #Convert amu to grams (1 amu = 1.66054e-24 g)

    total_mass = sum([atom_mass_gram[i]*st.session_state.num_atoms_arr[i] for i in range(np.size(atom_mass_gram))])
    volume = total_mass / density

    box_cm = volume**(1/3) #(cm)
    ##convert box lenght  in cm to armstrong
    box_arm = box_cm/1e-8
    return box_arm, volume


#def uniformPore(*, max_pore_size=st.session_state.max_pore_size_value):
def uniformPore(max_pore_size,box):
    '''
    This function obtains a uniform pore distribution
    The max_pore_size paramter can be changed to desired maximum pore size
    '''
    return np.random.rand() * (box*max_pore_size)


#def chiPore(*, max_pore_size = st.session_state.max_pore_size_value , df = 1,loc = 0, scale = .25, size=1):
def chiPore(max_pore_size, box, df = 1, size = 1):
    '''
    This function obtains a chi distributed pore distribution (https://en.wikipedia.org/wiki/Chi_distribution)
    The default parameter aguments are set to resemble a half-Gaussian
    The scale of .25 is set so that the maximum is around 1, to get the max allowed pore sized
    '''
    return np.random.chisquare(df, size)[0] * (box*max_pore_size)


#def betaPore(*, max_pore_size=st.session_state.max_pore_size_value , a = 4, b = 4):
def betaPore(max_pore_size, box, a = 4, b = 4):
    '''
    This function obatins a beta distribution (https://en.wikipedia.org/wiki/Beta_distribution)
    The parameters, a and b, set at 4 give a non-negative bell-shaped distribution, akin to a Gaussian Distribution
    '''
    return np.random.beta(a,b) * (box*max_pore_size)

#def poreCreator(*, porosity = st.session_state.porosity, num_pores=st.session_state.num_pores, epsilon=0.1, pore_dist_kind = 3):
def poreCreator(max_pore_size, box, porosity, num_pores,  pore_dist_kind = 3):

    '''
    This function creates the pores in a desired distribution and porosity
    
    porosity is the desired porosity, example; 0.5
    num_pores is the number of pores you want to sample
    epsilon ensures that we select a distribution that is less the desired porosity. We set this to 0.1
    pore_dist_kind specify what pore distribution is preffered example 3 (default) is the beta distribution
    '''
    epsilon=0.1
    poreRadii_list = []
    poreVolume_sum = 0
    i = 0           #iterator for our loop
    restarts = 0    #tells how many times this function is restarted when the distribution surpases the desired porsity
    porosity_threshold = porosity - epsilon


    while i < num_pores:
    
        if poreVolume_sum <= porosity_threshold:
            
            if pore_dist_kind == 1:
                name = "Beta Distribution"
                pore_radius = betaPore(max_pore_size,box)
                
            elif pore_dist_kind == 2:
                name = "Uniform Distribution"
                pore_radius = uniformPore(max_pore_size,box)
                
            else:
                name = "Chi Distribution"
                pore_radius = chiPore(max_pore_size,box)

            
            poreVolume_sum += pore_radius**3/box**3
            poreRadii_list.append(pore_radius)
            i += 1
        
              
            if poreVolume_sum > porosity:
                restarts += 1
                #print(f"Number of restarts = {restarts}", end = '\r')
                poreVolume_sum = 0
                i = 0
                poreRadii_list = []
            
            
            # Ensure the last pore sums to the desired porosity
            if i == num_pores-1:
                #print(f"Sampled up to {i} pores. Current Porosity is: {poreVolume_sum:.2f}")
                last_poreVolume = abs((porosity - poreVolume_sum)) * box**3
                #print(f"Added porosity fraction is {(last_poreVolume/box**3):.2f}")
                poreRadii_list.append(last_poreVolume**(1/3))
                poreVolume_sum += last_poreVolume/box**3
                #print(f"Added last pore. Final porosity is: {poreVolume_sum:.2f}")
                break
            
   
        if poreVolume_sum > porosity_threshold:
            restarts += 1
            poreVolume_sum = 0
            i = 0
            poreRadii_list = []
        
    #print()
    #print(f"                 *** Pore Radii List for {len(poreRadii_list)} Pores ***")
    #print(np.round(poreRadii_list,2))
    
    return name, poreRadii_list
##########################################################################################################


################################### MAIN ################################################################

col1, col2 = st.columns(2)


with col1:
    # Input parameters for pore distribution
    st.header('Pore Distribution')
    #pore_dist_kind = st.selectbox('Pore Distribution Kind', options=[1, 2, 3], format_func=lambda x: ["Uniform", "Chi", "Beta"][x-1])
    pore_dist_kind = st.radio('Select Pore Distribution', options=[1, 2, 3], format_func=lambda x: ["Beta", "Uniform", "Chi" ][x-1], horizontal=1)
    
    if st.button('Create Pores',key="createButton",on_click=disable, args=(False,)):
        st.session_state.box, volume = boxSize(st.session_state.density)
        st.session_state.name, st.session_state.poreRadii_list = poreCreator(st.session_state.max_pore_size, st.session_state.box, st.session_state.porosity, st.session_state.num_pores, pore_dist_kind)

       
        st.write(f"The box lenght is {st.session_state.box:.2f} \u212B")
        st.write(f"Pore Distribution: {st.session_state.name}")
        st.text(f"Pore Radii List:\n {np.round(st.session_state.poreRadii_list, 2)}")
        
        # Create distplot with custom bin_size

        st.session_state.df = pd.DataFrame(st.session_state.poreRadii_list, columns=["Pore Size"])
        
        st.session_state.fig_poreDistro = px.histogram(st.session_state.df, x = "Pore Size", nbins=int(st.session_state.num_pores/2))
        # Plot!
        st.plotly_chart(st.session_state.fig_poreDistro, use_container_width=True,theme="streamlit")

        
        
with col2:
    st.header("Generate Model")

    if st.button('Generate Model', key="generateButton",disabled=st.session_state.get("disabled", True)):
        with col1:
            st.write(f"The box lenght is {st.session_state.box:.2f} \u212B")
            st.write(f"Pore Distribution: {st.session_state.name}")
            st.text(f"Pore Radii List:\n {np.round(st.session_state.poreRadii_list, 2)}")
            # Create distplot with custom bin_size

            st.plotly_chart(st.session_state.fig_poreDistro, use_container_width=True ,theme="streamlit")


        ######################## Carbon foam constructor Algorithm starts here ####################################################

        pos = np.zeros([st.session_state.num_atoms+st.session_state.num_pores,3],float)   ## A list that takes in the position of the atoms and pores

        #st.write(f"Number of pores and atoms: {np.shape(pos)[0]} in {np.shape(pos)[1]} dimensions")
        
        box = st.session_state.box
        pore_overlap = st.session_state.pore_overlap
        poreRadii_list = st.session_state.poreRadii_list
        num_pores = st.session_state.num_pores
        num_atoms = st.session_state.num_atoms
        cutoff = st.session_state.cutoff

        ct = 0
        rnd = lambda i: i-round(i/box)*box  ## This makes rnd a function that takes i and perform i - round(i/box)*box
        vec_rnd = np.vectorize(rnd)         ## takes  a function and gives result in a callable vectorised function

        #print ("Now creating the center of the foams")

        while ct < num_pores:

            center =np.array([box*ran.random(),box*ran.random(),box*ran.random()])
            test = 0
            while test < ct:
                pore_distance = center-pos[test][:]
                pore_distance = vec_rnd(pore_distance)

                if sum(map(lambda i: i*i, pore_distance)) < (pore_overlap*poreRadii_list[test])**2:
                    center =np.array([box*ran.random(),box*ran.random(),box*ran.random()])
                    test = 0
                else:
                    test +=1 
           
            pos[ct][:] = center[:]
        
            ct += 1
    
      
        st.write("ALL PORES CREATED")

        ### Uncomment for debugging
        #print(pos[:num_pores])
        ###########################

        #print()
        with st.spinner("Creating carbon atoms. Takes time, please wait.)"):
            my_bar = st.progress(0, text="Progress Status.")

            while ct < num_atoms+num_pores:
                atoms =[box*ran.random(),box*ran.random(),box*ran.random()]
                test = 0
    
                ### Uncomment for debugging
                #print (atoms-pos[test][:])
                ###########################
    
                while test < ct:
                    atom_distance = atoms-pos[test][:]
                    atom_distance = vec_rnd(atom_distance)
        
                    # Make sure atoms are not close to the pore center
                    if test < num_pores-1 and sum(map(lambda i: i*i, atom_distance)) < poreRadii_list[test]**2:
                        atoms =np.array([box*ran.random(),box*ran.random(),box*ran.random()])
            
                        ### Uncomment for debugging##################
                        #print ("OOPS!!! TOO CLOSE to a foam center")
                        #############################################
                        test = 0 
            
                    ### Position the atoms if atoms are not close to center    
                    elif sum(map(lambda i: i*i, atom_distance)) < cutoff**2 and test >= num_pores:
                        atoms =np.array([box*ran.random(),box*ran.random(),box*ran.random()])
                        
                        ### Uncomment for debugging##################
                        #print ("OOPS!!! TOO CLOSE to a foam center")
                        #############################################
            
                        test = 0
                    else:
                        test += 1
            
                pos[ct][:] = atoms
    
    
                ct +=1
                my_bar.progress(ct/(num_atoms+num_pores), text="Progress Status.")         
                #print (f"Placing Atom number {ct-num_pores} of {num_atoms}", end='\r')
         ################################################################################################

        my_bar.empty()
        with st.spinner("Preparing POSCAR file for download. Please wait..."):
            time.sleep(5)

    
        atom_position = pos/st.session_state.box

        #Convert the NumPy array to a list of lists
        list_of_lists = atom_position.tolist()

        #Convert each sublist to a string with elements separated by spaces
        lines = [" ".join(map(str, sublist)) for sublist in list_of_lists]

        #Join the resulting lines with newline characters
        final_output = "\n".join(lines)

        st.session_state.txt1 = st.text_area("POSCAR File", f"{stringNumAtoms} {stringDensity} {stringNoFoam} {stringFoamOverlap}\n\
{st.session_state.box:10.6f}\n\
{1.0:2.6f} {0.0:2.6f} {0.0:2.6f}\n\
{0.0:2.6f} {1.0:2.6f} {0.0:2.6f}\n\
{0.0:2.6f} {0.0:2.6f} {1.0:2.6f}\n\
O     C \n\
{st.session_state.num_pores}   {st.session_state.num_atoms} \n\
Direct\n{final_output}")
                                               


        st.success('Done!')
        st.info(f'IMPORTANT! The first {st.session_state.num_pores} coordinates in the POSCAR file below designate the pores as oxygen (O) atoms.', icon="ℹ️")

        if st.download_button(label="Download POSCAR",data=st.session_state.txt1, file_name=st.session_state.atoms_vasp,on_click=disable, args=(True,)):
            st.write("Download Complete.")
