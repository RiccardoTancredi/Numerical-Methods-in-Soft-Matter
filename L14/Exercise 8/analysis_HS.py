import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import matplotlib.gridspec as gridspec
import os

plt.rcParams.update({'font.size': 15})  
folder = 'res/'

# plt.plot(prey, predators, label=r'$\mathcal{C}(0) = ' + f'{tuple(data[0])}' + r'$', 
#          color='cornflowerblue')
# plt.plot(prey_1, predators_1, label=r'$\mathcal{C}(0) = ' + f'{tuple(data_1[0])}' + r'$', 
#          color='teal')
# plt.ylim((0, 1500))
# plt.xlabel('Prey')
# plt.ylabel('Predators')
# plt.title('Lotka - Volterra')
# plt.legend()
# plt.grid(ls='--')
# plt.tight_layout()
# # plt.savefig(f'img/prey_predators_dynamics.png', dpi=300)
# plt.show()


#######################
# MC initial position #
#######################

initial = False
if initial:
    read_file = np.loadtxt('SOURCE/param.dat', dtype=str)
    var_names, data = [], []
    for ar in read_file:
        for el in ar:
            try:
                data.append(int(el))
            except:
                if '.dat' not in el:
                    var_names.append(el)
                else:
                    if 'velocity' not in el:
                        start_file_name = str(el)
                    else:
                        velocity_file_name = str(el)                
                    
    var_names = np.array(var_names)
    var_names = var_names[var_names != '=']
    N, N_strep, box_x, box_y, box_z, T, sigma = data[:7]
    iteration = data[-1]
    # k_B = 1.38e-23
    # m = k_B*T*sigma**2

    file_path_x = f'SOURCE/{start_file_name}'
    file_path_v = f'SOURCE/{velocity_file_name}'
    if not os.path.exists(file_path_x) or not os.path.exists(file_path_v):
        # Generate initial (x, y, z) positions of the N particles
        init_x = np.random.uniform(low=0, high=box_x, size=(N, 3)) 
        init_v = np.random.normal(loc=0, scale=sigma, size=(N, 3))
        np.savetxt(file_path_x, init_x)
        np.savetxt(file_path_v, init_v)
        print('Files created!')
    else:
        init_x = np.loadtxt(file_path_x)
        init_v = np.loadtxt(file_path_v)
        print('Files already exist.')


show_plots = False
####################
# 3d visualization #
####################
    
if show_plots:

    with open('SOURCE/res/restartpoint.dat', 'r') as file:
        lines = file.readlines()
    data = []
    for line in lines:
        if line.strip():  # Skip empty lines
            x, y, z = map(float, line.split())
            data.append([x, y, z])

    data = np.array(data)
    energy = np.loadtxt('SOURCE/res/energy.dat')

    # energy plot
    plt.plot(energy, ls='-.', label=r'$E$', color='royalblue')
    plt.grid(ls='--')
    plt.legend()
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\log(E)$')
    plt.yscale('log')
    plt.tight_layout()
    plt.show()

    # Create 3d plot:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    particles, = ax.plot([], [], [], 'bo', markersize=6)


    def update(frame):
        ax.cla()  # Clear previous frame
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')
        ax.set_zlabel(r'$z$')
        ax.set_title(f'Frame {frame}') # , E = {energy[frame]}

        ax.set_xlim([0, box_x])
        ax.set_ylim([0, box_y])
        ax.set_zlim([0, box_z])

        # Plot the particles for the current frame
        particles, = ax.plot(data[frame * N:(frame + 1) * N][:, 0],
                            data[frame * N:(frame + 1) * N][:, 1],
                            data[frame * N:(frame + 1) * N][:, 2], 'bo', markersize=6)

        # if frame == num_frames - 1:
        #     animation.event_source.stop()

    # Set the number of frames based on the number of frames in your data
    num_frames = len(data) // N

    # Create the animation
    animation = FuncAnimation(fig, update, frames=num_frames, interval=100, blit=False)

    plt.show()


    ###################
    # 2d scatter plot #
    ###################

    fig, ax = plt.subplots()
    scatter = ax.scatter([], [], c=[], cmap='viridis', marker='o')

    # Set the axis labels and title
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_title('Heatmap Animated Plot')

    # Get the initial color limits
    initial_clim = scatter.get_clim()*(box_z)

    # Update function for animation
    def update(frame):
        ax.set_title(f'Frame {frame}')
        x = data[frame * N:(frame + 1) * N][:, 0]
        y = data[frame * N:(frame + 1) * N][:, 1]
        z = data[frame * N:(frame + 1) * N][:, 2]

        ax.set_xlim([0, box_x])
        ax.set_ylim([0, box_y])

        # Update the scatter plot
        scatter.set_offsets(np.column_stack((x, y)))
        scatter.set_array(z/box_z)


    # Create the animation
    animation = FuncAnimation(fig, update, frames=num_frames, interval=500, blit=False)

    # Set colorbar
    cbar = plt.colorbar(scatter)
    cbar.set_label(r'$z$')

    # Show the plot
    plt.show()


for keyword in ['random', 'cube']:   # 
    delta_x = ['0.01', '0.05', '0.3', '0.5', '1']
    files = ['acceptance.dat', 'energy.dat', 'restartpoint.dat']
    folder = f'res/{keyword}/'
    save_fig = f'../img/{keyword}_'
    colors = np.array(["#ddfff7", "#93e1d8", "#ffa69e", "#aa4465", "#460b2e"])
    colors = colors[::-1]
    with open('param.dat') as file:
        lines = file.readlines()
    N = int(lines[0].split('=')[-1])
    rho = np.array([0.05, 0.3, 0.5, 1][::-1])
    volum = np.round((N/rho)**(1/3)).astype(int)
    gs = gridspec.GridSpec(2, 3, width_ratios=[1, 1, 1], height_ratios=[1, 1])

    acceptance_plot = True

    if acceptance_plot:

        for k, V in enumerate(volum):
            print(f'Buffering rho = {rho[k]} ...')
            fig, axs = plt.subplots(2, 3, figsize=(12, 8), gridspec_kw={'width_ratios': [1, 1, 1], 'height_ratios': [1, 1]})
            min_y, max_y = 0, 0
            for i, d in enumerate(delta_x):
                row, col = i // 3, i % 3 
                data = []

                for it in range(1, 11):
                    file_name = f'{folder}{V}_{d}/{it}/V_{V}_'
                    # here we deal only with acceptance
                    file_name += files[0]
                    data.append(np.loadtxt(file_name))
                
                N_steps = len(data[0])
                data = np.array(data)/N_steps/N # -> rate
                means = data.mean(axis=0)
                # print(data.shape, means.shape)
                yerr = means.std(axis=0)/np.sqrt(10)  # std associated to the mean

                axs[row, col].errorbar(x=list(range(len(means))), y=means, 
                                    yerr=yerr, color=colors[i], fmt='.',
                                    label=r'$d_{max}=$'+f'{float(d)}')
                min_y = min(min_y, np.min(means-yerr))
                max_y = max(max_y, np.max(means+yerr))

                # axs[row, col].scatter(x=list(range(len(means))), y=means, 
                #                       color=colors[i], marker='.',
                #                       label=r'$d_{max}=$'+f'{float(d)}')
                # axs[row, col].legend(loc='upper left')
                axs[row, col].grid(ls='--', alpha=0.5)
                axs[row, col].set_xscale('log')
                # axs[row, col].set_yscale('log')
                axs[row, col].set_xlabel(r'$\log(t)$')
                if col == 0:
                    axs[row, col].set_ylabel('Acceptance rate')

            fig.delaxes(axs[1, 2])
            
            handles, labels = [], []
            for ax in axs.flat:
                handles_, labels_ = ax.get_legend_handles_labels()
                handles.extend(handles_)
                labels.extend(labels_)

            leg = fig.legend(handles, labels, loc='center', bbox_to_anchor=(0.85, 0.3), borderaxespad=0.)

            # Set the same y-axis limits for all subplots
            for ax in axs.flat:
                ax.set_ylim([float(min_y), float(max_y)])

            # plt.suptitle(r'$\rho = ' + f'{rho[k]}' + r'$')
            plt.tight_layout()
            plt.savefig(f'{save_fig}rho_{rho[k]}.png', dpi=300, bbox_extra_artists=(leg,), bbox_inches='tight')
            plt.close()
            # plt.show()

    energy_plot = True

    if energy_plot:
        
        for k, V in enumerate(volum):
            print(f'Buffering energy - rho = {rho[k]} ...')
            fig, axs = plt.subplots(2, 3, figsize=(12, 8), gridspec_kw={'width_ratios': [1, 1, 1], 'height_ratios': [1, 1]})
            
            min_y, max_y = 0, 0
            for i, d in enumerate(delta_x):
                row, col = i // 3, i % 3 
                data = []

                for it in range(1, 11):
                    file_name = f'{folder}{V}_{d}/{it}/V_{V}_'
                    # here we deal only with acceptance
                    file_name += files[1]
                    data.append(np.loadtxt(file_name))
                
                data = np.array(data)
                means = data.mean(axis=0)
                # print(data.shape, means.shape)
                yerr = means.std(axis=0)/np.sqrt(10)  # std associated to the mean

                # np.savetxt(f'{folder}{V}_means_en.dat', means)
                # np.savetxt(f'{folder}{V}_yerr_en.dat', yerr)
                # else:
                #     means = np.loadtxt(f'{folder}{V}_means_en.dat')
                #     yerr = np.loadtxt(f'{folder}{V}_yerr_en.dat')

                axs[row, col].errorbar(x=list(range(len(means))), y=means, 
                                       yerr=yerr, color=colors[i], fmt='.',
                                       label=r'$d_{max}=$'+f'{float(d)}')
                # axs[row, col].scatter(x=list(range(len(means))), y=means, 
                #                     color=colors[i], marker='.',
                #                     label=r'$d_{max}=$'+f'{float(d)}')
                
                min_y = min(min_y, np.min(means-yerr))
                max_y = max(max_y, np.max(means+yerr))

                # axs[row, col].legend(loc='upper left')
                axs[row, col].grid(ls='--', alpha=0.5)
                axs[row, col].set_xscale('log')
                # axs[row, col].set_yscale('log')
                axs[row, col].set_xlabel(r'$\log(t)$')
                if col == 0:
                    axs[row, col].set_ylabel('Energy')

            fig.delaxes(axs[1, 2])
            
            handles, labels = [], []
            for ax in axs.flat:
                handles_, labels_ = ax.get_legend_handles_labels()
                handles.extend(handles_)
                labels.extend(labels_)

            leg = fig.legend(handles, labels, loc='center', bbox_to_anchor=(0.85, 0.3), borderaxespad=0.)

            # Set the same y-axis limits for all subplots
            for ax in axs.flat:
                ax.set_ylim([-100, float(max_y)])

            # plt.suptitle(r'$\rho = ' + f'{rho[k]}' + r'$')
            plt.tight_layout()
            plt.savefig(f'{save_fig}energy_rho_{rho[k]}.png', dpi=300, bbox_extra_artists=(leg,), bbox_inches='tight')
            plt.close()
            # plt.show()
