import numpy as np
from os import path
import pickle
from ev_model import persistency

def autocorrFFT(x):
    N=len(x)
    F = np.fft.fft(x, n=2*N)  #2*N because of zero-padding
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD)
    res= (res[:N]).real   #now we have the autocorrelation in convention B
    n=N*np.ones(N)-np.arange(0,N) #divide res(m) by (N-m)
    return res/n #this is the autocorrelation in convention A

def msd_per_radius(displacements):
    displacement_per_radius = {}
    for k,v in displacements.items():
        if displacements[k].radius_um not in displacement_per_radius:
            displacement_per_radius[displacements[k].radius_um] = []
        displacement_per_radius[displacements[k].radius_um].append(displacements[k].msd_fft())
    for k, v in displacement_per_radius.items():
        displacement_per_radius[k] = np.array(displacement_per_radius[k]).mean()
    return displacement_per_radius # pd.DataFrame.from_dict(displacement_per_radius)

def obtain_evs_per_radius(displacements):
    evs_per_radius = {}
    for i in range(len(displacements)):
        if displacements[i].radius_um in evs_per_radius:
            evs_per_radius[displacements[i].radius_um].append(i)
        else:
            evs_per_radius[displacements[i].radius_um] = [i]
    return evs_per_radius

def calculate_diffusion_rate(ev, dt, tao):
    X = np.array([ev.x, ev.y])
    # compute the instantaneous velocity at each time point
    inst_velocities = np.diff(X) / (dt * tao)
    # compute the variance
    velocity_variance = np.var(inst_velocities, axis=1, ddof=1)
    # velocity_variance
    mean_velocity_variance = np.mean(velocity_variance)
    #mean_velocity_variance
    estimated_D = mean_velocity_variance / 2 * (dt * tao)
    return estimated_D

def compute_diffusion_rate_per_radius(displacements, evs_per_radius, dt, tao):
    # sort the keys
    sorted_ev_radius = [k for k in evs_per_radius.keys()]
    sorted_ev_radius.sort()

    diffusion_rates = {'radius_in_um':[], 'diffusion_rate_mean':[], 'diffusion_rate_sd':[], 'analytical':[]}#np.zeros(len(sorted_ev_radius))
    for radius in sorted_ev_radius:
        diffusions = []
        for ev_id in evs_per_radius[radius]:
            # compute the speed and accumulate the values
            diffusion_rate = calculate_diffusion_rate(displacements[ev_id], dt, tao)
            diffusions.append(diffusion_rate)
        # now we compute the mean and sd
        diffusion_rates['radius_in_um'].append(radius)
        diffusion_rates['diffusion_rate_mean'].append(np.mean(diffusions))
        diffusion_rates['diffusion_rate_sd'].append(np.std(diffusions))
        diffusion_rates['analytical'].append((4.05E-3/(0.018879715 * radius)))  # E-9
    return diffusion_rates

class EV_positions:
    def __init__(self, _id, radius_um, x, y):
        self._id = _id
        self.radius_um = radius_um
        self.x = [x]
        self.y = [y]
        self.iteration = []

    def compute_displacement(self):
        self.dx = self.x[1:] - self.x[:-1]
        self.dy = self.y[1:] - self.y[:-1]
        return [self.dx, self.dy]

    def compute_msd(self):
        acum = 0
        n = len(self.x)
        for i in range(1, n):
            distance = np.sqrt((self.x[i] - self.x[i-1])**2 + (self.y[i] - self.y[i-1])**2)
            acum += distance
        self.msd = acum / n
        return self.msd

    def compute_msd_straight_forward(self):
        """
        https://stackoverflow.com/a/34222273
        """
        rows = np.array((self.x, self.y)).T
        shifts = np.arange(len(rows))
        msds = np.zeros(shifts.size)

        for i, shift in enumerate(shifts):
            diffs = rows[:-shift if shift else None] - rows[shift:]
            sqdist = np.square(diffs).sum(axis=1)
            msds[i] = sqdist.mean()
        return msds

    def msd_fft(self):
        rows = np.array((self.x, self.y)).T
        N=len(rows)
        D=np.square(rows).sum(axis=1)
        D=np.append(D,0)
        S2=sum([autocorrFFT(rows[:, i]) for i in range(rows.shape[1])])
        Q=2*D.sum()
        S1=np.zeros(N)
        for m in range(N):
            Q=Q-D[m-1]-D[N-m]
            S1[m]=Q/(N-m)
        return S1-2*S2

def obtain_displacements_at_intervals(base_path, pts0, top_limit, tao, overwrite=False):
    """
    Obtains the RMS displacement per EV.
    Checks if the information is already available in pickled format and loads it.
    Otherwise, it computes the values and saves them in pickled format for future re use.
    This function only considers those EVs existing in the file provided in pts0 which holds the first state of the simulation.

    Inputs:
    pts0 - the XML file with the first state of the model we should read
    top_limit - the desired iteration number to fetch the last state of the system from from the saved XML files.
    tao - the lag or stride to use for computing the RMS value

    Outputs:
    displacements - A list of EV objects
    """
    displacements = {}

    target_displacements_file = path.join(base_path, f"{(pts0[:-4] if pts0.endswith('.xml') else pts0)}_displ_at_{tao}.pickle")
    if path.exists(target_displacements_file) and not overwrite:
        with open(target_displacements_file, 'rb') as target_file:
            print(f'Loading values from pickled file {target_displacements_file}')
            return pickle.load(target_file)
    else:
        # Read the iterations, compute the values, and save to file
        print("Computing values...")
        for iteration in range(0, top_limit, tao):
            file_name = path.join(base_path, (pts0 if iteration == 0 else f'{iteration}.xml'))
            if path.isfile(file_name):
                #read_iteration_data_to_columns(basedir, states_file, dt = -1, prefix='')
                itno, environment, secretory, ciliary, evs = persistency.read_xml(base_path, file_name)
                parameters = ['id', 'x', 'y', 'radius_um', 'age']
                evs_df = persistency.agent_data_to_data_frame(evs, parameters)

                # Only those EVs existing in 0.xml should be considered
                for idx in range(evs_df.shape[0]):
                    evid = evs_df.iloc[idx].id
                    if evid in displacements:
                        displacements[evid].x.append(evs_df.iloc[idx].x)
                        displacements[evid].y.append(evs_df.iloc[idx].y)
                        displacements[evid].iteration.append(iteration)
                    else:
                        if iteration == 0:
                            displacements[evid] = EV_positions(evs_df.iloc[iteration].id, 
                                            evs_df.iloc[iteration].radius_um, 
                                            evs_df.iloc[iteration].x, evs_df.iloc[iteration].y)
            else:
                raise ValueError(f'It\'s NOT a file! {file_name}')
        # save the computed values
        print('Saving the computed values as:', target_displacements_file)
        with open(target_displacements_file, 'wb') as target_output:
            pickle.dump(displacements, target_output)
            print("Done!")
    return displacements
