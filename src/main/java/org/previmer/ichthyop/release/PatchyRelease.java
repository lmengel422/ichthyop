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

import org.previmer.ichthyop.TypeZone;
import org.previmer.ichthyop.Zone;
import org.previmer.ichthyop.particle.IParticle;
import org.previmer.ichthyop.event.ReleaseEvent;
import org.previmer.ichthyop.io.ZoneTracker;
import org.previmer.ichthyop.particle.ParticleFactory;

/**
 *
 * @author pverley
 */
public class PatchyRelease extends AbstractRelease {

    private int nb_patches;
    private int nb_agregated;
    private double radius_patch, thickness_patch;
    private int nbReleaseZones;
    private boolean is3D;
    private boolean perZone;
    private static final double ONE_DEG_LATITUDE_IN_METER = 111138.d;

    @Override
    public void loadParameters() throws Exception {

        /* retrieve patches parameters */
        nb_patches = Integer.valueOf(getParameter("number_patches"));
        nb_agregated = Integer.valueOf(getParameter("number_agregated"));
        radius_patch = Float.valueOf(getParameter("radius_patch"));
        thickness_patch = Float.valueOf(getParameter("thickness_patch"));
        // new feature, the parameter does not exist in ichthyop <= 3.2
        try {
            perZone = Boolean.valueOf(getParameter("per_zone"));
        } catch (Exception ex) {
            perZone = false;
        }

        /* Check whether 2D or 3D simulation */
        is3D = getSimulationManager().getDataset().is3D();

        /* Load release zones*/
        getSimulationManager().getZoneManager().loadZonesFromFile(getParameter("zone_file"), TypeZone.RELEASE);
        nbReleaseZones = (null != getSimulationManager().getZoneManager().getZones(TypeZone.RELEASE))
                ? getSimulationManager().getZoneManager().getZones(TypeZone.RELEASE).size()
                : 0;
        getSimulationManager().getOutputManager().addPredefinedTracker(ZoneTracker.class);
    }

    @Override
    public int release(ReleaseEvent event) throws Exception {
        if (perZone) {
            return releasePerZone();
        } else {
            return releaseUniform();
        }
    }

    private int releasePerZone() {

        double xmin, xmax, ymin, ymax;
        double upDepth, lowDepth;

        int index = Math.max(getSimulationManager().getSimulation().getPopulation().size(), 0);
        for (int i_zone = 0; i_zone < nbReleaseZones; i_zone++) {
            Zone zone = getSimulationManager().getZoneManager().getZones(TypeZone.RELEASE).get(i_zone);
            xmin = zone.getXmin();
            xmax = zone.getXmax();
            ymin = zone.getYmin();
            ymax = zone.getYmax();
            upDepth = zone.getUpperDepth();
            lowDepth = zone.getLowerDepth();
            for (int p = 0; p < nb_patches; p++) {
                int DROP_MAX = 2000;
                IParticle particle = null;
                int counter = 0;
                while (null == particle) {
                    if (counter++ > DROP_MAX) {
                        throw new NullPointerException("{Patchy release} Unable to release particle. Check out the zone definitions.");
                    }
                    double x = xmin + Math.random() * (xmax - xmin);
                    double y = ymin + Math.random() * (ymax - ymin);
                    double depth = Double.NaN;
                    if (is3D) {
                        depth = -1.d * (upDepth + Math.random() * (lowDepth - upDepth));
                    }
                    particle = ParticleFactory.createZoneParticle(index, x, y, depth);
                }
                getSimulationManager().getSimulation().getPopulation().add(particle);
                index++;
                for (int f = 0; f < nb_agregated - 1; f++) {
                    IParticle particlePatch = null;
                    counter = 0;
                    while (null == particlePatch) {

                        if (counter++ > DROP_MAX) {
                            throw new NullPointerException("{Patchy release} Unable to release particle. Check out the patchy release definition.");
                        }
                        double lat = particle.getLat() + radius_patch * (Math.random() - 0.5d) / ONE_DEG_LATITUDE_IN_METER;
                        double one_deg_longitude_meter = ONE_DEG_LATITUDE_IN_METER * Math.cos(Math.PI * particle.getLat() / 180.d);
                        double lon = particle.getLon() + radius_patch * (Math.random() - 0.5d) / one_deg_longitude_meter;
                        double depth = Double.NaN;
                        if (is3D) {
                            depth = particle.getDepth() + thickness_patch * (Math.random() - 0.5d);
                        }
                        particlePatch = ParticleFactory.createGeoParticle(index, lon, lat, depth);
                    }
                    getSimulationManager().getSimulation().getPopulation().add(particlePatch);
                    index++;
                }
            }
        }
        return index;

    }

    private int releaseUniform() {

        double xmin, xmax, ymin, ymax;
        double upDepth = Double.MAX_VALUE, lowDepth = 0.d;
        /**
         * Reduces the release area function of the user-defined zones
         */
        xmin = Double.MAX_VALUE;
        ymin = Double.MAX_VALUE;
        xmax = 0.d;
        ymax = 0.d;
        for (int i_zone = 0; i_zone < nbReleaseZones; i_zone++) {
            Zone zone = getSimulationManager().getZoneManager().getZones(TypeZone.RELEASE).get(i_zone);
            xmin = Math.min(xmin, zone.getXmin());
            xmax = Math.max(xmax, zone.getXmax());
            ymin = Math.min(ymin, zone.getYmin());
            ymax = Math.max(ymax, zone.getYmax());
            if (is3D) {
                upDepth = Math.min(upDepth, zone.getUpperDepth());
                lowDepth = Math.max(lowDepth, zone.getLowerDepth());
            } else {
                upDepth = lowDepth = Double.NaN;
            }
        }

        int index = Math.max(getSimulationManager().getSimulation().getPopulation().size(), 0);
        for (int p = 0; p < nb_patches; p++) {
            /**
             * Instantiate a new Particle
             */
            int DROP_MAX = 2000;
            IParticle particle = null;
            int counter = 0;
            while (null == particle) {
                if (counter++ > DROP_MAX) {
                    throw new NullPointerException("{Patchy release} Unable to release particle. Check out the zone definitions.");
                }
                double x = xmin + Math.random() * (xmax - xmin);
                double y = ymin + Math.random() * (ymax - ymin);
                double depth = Double.NaN;
                if (is3D) {
                    depth = -1.d * (upDepth + Math.random() * (lowDepth - upDepth));
                }
                particle = ParticleFactory.createZoneParticle(index, x, y, depth);
            }
            getSimulationManager().getSimulation().getPopulation().add(particle);
            index++;
            for (int f = 0; f < nb_agregated - 1; f++) {

                IParticle particlePatch = null;
                counter = 0;
                while (null == particlePatch) {

                    if (counter++ > DROP_MAX) {
                        throw new NullPointerException("{Patchy release} Unable to release particle. Check out the patchy release definition.");
                    }
                    double lat = particle.getLat() + radius_patch * (Math.random() - 0.5d) / ONE_DEG_LATITUDE_IN_METER;
                    double one_deg_longitude_meter = ONE_DEG_LATITUDE_IN_METER * Math.cos(Math.PI * particle.getLat() / 180.d);
                    double lon = particle.getLon() + radius_patch * (Math.random() - 0.5d) / one_deg_longitude_meter;
                    double depth = Double.NaN;
                    if (is3D) {
                        depth = particle.getDepth() + thickness_patch * (Math.random() - 0.5d);
                    }
                    particlePatch = ParticleFactory.createGeoParticle(index, lon, lat, depth);
                }
                getSimulationManager().getSimulation().getPopulation().add(particlePatch);
                index++;
            }
        }
        return index;
    }

    @Override
    public int getNbParticles() {
        int multiplier = perZone ? getSimulationManager().getZoneManager().getZones(TypeZone.RELEASE).size() : 1;
        return multiplier * nb_patches * nb_agregated;
    }
}
