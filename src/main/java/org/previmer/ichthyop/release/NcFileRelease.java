/*
 *
 * ICHTHYOP, a Lagrangian tool for simulating ichthyoplankton dynamics
 * http://www.ichthyop.org
 *
 * Copyright (C) IRD (Institut de Recherce pour le Developpement) 2006-2020
 * http://www.ird.fr
 *
 * Main developper: Philippe VERLEY (philippe.verley@ird.fr), Nicolas Barrier (nicolas.barrier@ird.fr)
 * Contributors (alphabetically sorted):
 * Gwendoline ANDRES, Sylvain BONHOMMEAU, Bruno BLANKE, Timothee BROCHIER,
 * Christophe HOURDIN, Mariem JELASSI, David KAPLAN, Fabrice LECORNU,
 * Christophe LETT, Christian MULLON, Carolina PARADA, Pierrick PENVEN,
 * Stephane POUS, Nathan PUTMAN.
 *
 * Ichthyop is a free Java tool designed to study the effects of physical and
 * biological factors on ichthyoplankton dynamics. It incorporates the most
 * important processes involved in fish early life: spawning, movement, growth,
 * mortality and recruitment. The tool uses as input time series of velocity,
 * temperature and salinity fields archived from oceanic models such as NEMO,
 * ROMS, MARS or SYMPHONIE. It runs with a user-friendly graphic interface and
 * generates output files that can be post-processed easily using graphic and
 * statistical software.
 *
 * To cite Ichthyop, please refer to Lett et al. 2008
 * A Lagrangian Tool for Modelling Ichthyoplankton Dynamics
 * Environmental Modelling & Software 23, no. 9 (September 2008) 1210-1214
 * doi:10.1016/j.envsoft.2008.02.005
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 3 of the License). For a full
 * description, see the LICENSE file.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

package org.previmer.ichthyop.release;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.logging.Level;
import org.previmer.ichthyop.particle.IParticle;
import org.previmer.ichthyop.event.ReleaseEvent;
import org.previmer.ichthyop.io.IOTools;
import org.previmer.ichthyop.particle.ParticleFactory;
import org.previmer.ichthyop.particle.ParticleMortality;
import ucar.ma2.ArrayDouble;
import ucar.ma2.ArrayFloat;
import ucar.ma2.ArrayInt;
import ucar.nc2.NetcdfFile;
import ucar.nc2.dataset.NetcdfDatasets;

/**
 *
 * @author pverley
 */
public class NcFileRelease extends AbstractRelease {

    private String filename;

    @Override
    public void loadParameters() throws IOException {

        /* make sure file exits and can be read */
        filename = IOTools.resolveFile(getParameter("ncfile"));

        File file = new File(filename);
        if (!file.isFile()) {
            throw new FileNotFoundException("NetCDF file " + filename + " not found.");
        }
        if (!file.canRead()) {
            throw new IOException("NetCDF file " + file + " cannot be read");
        }
    }

    @Override
    public int release(ReleaseEvent event) throws Exception {

        double time = event.getSource().getTime();
        int index = Math.max(getSimulationManager().getSimulation().getPopulation().size(), 0);

        NetcdfFile nc = NetcdfDatasets.openFile(filename, null);
        ArrayDouble.D1 timeArr = (ArrayDouble.D1) nc.findVariable("time").read();
        int rank = 0;
        int length = timeArr.getShape()[0];

        /** Find current rank */
        double time_rank0 = timeArr.get(0), time_rank1;
        while (rank < length) {
            rank++;
            time_rank1 = timeArr.get(rank);
            if ((time >= time_rank0 && time < time_rank1)
                    || (time <= time_rank0 && time > time_rank1)) {
                break;
            }
            if (time == time_rank1) {
                rank++;
                time_rank0 = time_rank1;
                break;
            }
            time_rank0 = time_rank1;
        }
        rank--;

        double x = 0.f;
        if (time != time_rank0) {
            x = (time - timeArr.get(rank))
                    / (timeArr.get(rank + 1) - timeArr.get(rank));
        }

        ArrayFloat.D2 lonArr = (ArrayFloat.D2) nc.findVariable("lon").read();
        ArrayFloat.D2 latArr = (ArrayFloat.D2) nc.findVariable("lat").read();
        ArrayInt.D2 mortalityArr = (ArrayInt.D2) nc.findVariable("mortality").read();
        boolean bln3D = getSimulationManager().getDataset().is3D();
        ArrayFloat.D2 depthArr = null;
        if (bln3D) {
            depthArr = (ArrayFloat.D2) nc.findVariable("depth").read();
        }
        double lon, lat, depth = Double.NaN;
        IParticle particle;
        int nb_particles = lonArr.getShape()[1];

        for (int i = 0; i < nb_particles; i++) {

            if (x == 0.f) {
                lon = lonArr.get(rank, i);
                lat = latArr.get(rank, i);
                if (bln3D) {
                    depth = depthArr.get(rank, i);
                }
            } else {
                lon = x * lonArr.get(rank + 1, i)
                        + (1 - x) * lonArr.get(rank, i);
                lat = x * latArr.get(rank + 1, i)
                        + (1 - x) * latArr.get(rank, i);
                if (bln3D) {
                    depth = x * depthArr.get(rank + 1, i)
                            + (1 - x) * depthArr.get(rank, i);
                }
            }

            particle = ParticleFactory.createGeoParticle(index, lon, lat, depth, ParticleMortality.getMortality(mortalityArr.get(rank, i)));
            getSimulationManager().getSimulation().getPopulation().add(particle);
            index++;
        }

        nc.close();
        return index;
    }

    @Override
    public int getNbParticles() {

        try {
            NetcdfFile nc = NetcdfDatasets.openDataset(filename);
            int nParticle = nc.findDimension("drifter").getLength();
            nc.close();
            return nParticle;
        } catch (IOException ex) {
            getLogger().log(Level.SEVERE, null, ex);
            System.exit(1);
        }
        return -1;
    }
}
