import pandas as pd
import matplotlib.pyplot as plt
from math import sqrt, pi, e
import numpy as np

# Displays all columns of the Pandas DataFrame
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

# looks up spacial wind tunnel conditions
def look_up_flow_velocities(x, y, z):
    """returns x, y, and z wind tunnel flow velocities in m/s given the position of the droplet"""
    flow_x_velocity = 60.0
    flow_y_velocity = 0.0
    flow_z_velocity = 0.0
    return flow_x_velocity, flow_y_velocity, flow_z_velocity

def look_up_pressure(x, y, z):
    """returns pressure at a given position in Pa"""
    pressure = 101325
    return pressure


def look_up_air_temperature(x, y, z):
    """looks up air temperature in °C at a given position"""
    air_temperature = -13
    return air_temperature

def look_up_relative_humidity(x, y, z):
    """looks up relative humidity as a percent at a given position"""
    relative_humidity = 0.95
    return relative_humidity

# Looks up initial droplet conditions (could be turned into a function parameter)
def look_up_initial_droplet_velocity():
    """returns initial droplet velocity in m/s"""
    droplet_x_velocity = 51.0
    droplet_y_velocity = 0.0
    droplet_z_velocity = 0.0
    return droplet_x_velocity, droplet_y_velocity, droplet_z_velocity

def look_up_droplet_temperature():
    """returns droplet temperature in C°"""
    droplet_temperature = -3.5
    return droplet_temperature

# calculations
def calculate_droplet_volume(diameter):
    """Calculates the droplet volume in m^3 assuming a prefect sphere"""
    return pi * diameter**3 / 6

def calculate_film_temperature(air_temperature, T_surface):
    """calculates film temperature"""
    return (air_temperature + T_surface) / 2


def calculate_relative_velocity(flow_velocity, droplet_velocity ):
    """Calculates the relative velocity between the droplet velocity and the flow velocity"""
    return flow_velocity - droplet_velocity

def calculate_relative_velocity_magnitude(flow_x_velocity, flow_y_velocity, flow_z_velocity):
    """Calculates the magnitude of the relative velocity between the flow velocity and the droplet velocity"""
    v_rel = sqrt(flow_x_velocity**2 + flow_y_velocity**2 + flow_z_velocity**2)
    return v_rel

def calculate_water_density(T_mean):
    """Calculates the water density in kg/m3 using a curve fit equation"""
    water_density = (
              999.83952
              + 16.945176 * T_mean
              - 7.9870401e-3 * T_mean** 2
              - 46.170461e-6 * T_mean ** 3
              + 105.56302e-9 * T_mean** 4
              - 280.54253e-12 * T_mean ** 5
          ) / (1 + 16.879850e-3 * T_mean)
    return water_density

def calculate_air_density(pressure, specific_gas_constant_of_air, air_temperature):
    """ Calculates the air density in kg/m3 using the ideal gas law"""
    air_density = pressure/(specific_gas_constant_of_air * (air_temperature + 273.15))
    return air_density

def calculate_dynamic_viscosity(T_film):
    """Takes film temperature in °C and returns dynamic viscosity in kg/m3 using a curve fit equation"""
    dynamic_viscosity = (1.716 * 10 ** -5) * (((T_film +273.15)/273.15)**(3/2)) * ((273.15 + 111)/(T_film + 273.15 + 111))
    return dynamic_viscosity

def calculate_reynolds_number(air_density, v_rel, d, dynamic_viscosity):
    """Calculates Reynolds number"""
    reynolds_number = (air_density * v_rel * d) / dynamic_viscosity
    return reynolds_number

def calculate_drag_coefficient(reynolds_number):
    """Calculates drag coefficient using the Shiller-Naumann correlation"""
    if reynolds_number < 1000:
        drag_coefficient = (24/reynolds_number) * (1 + 0.15 * reynolds_number ** 0.687)
    else:
        drag_coefficient = max((24/reynolds_number) * (1 + 0.15 * reynolds_number ** 0.687), 0.44)
    return drag_coefficient

def calculate_projected_area(d,):
    """Calculates projected area assuming the droplet is a perfect sphere"""
    projected_area = pi * (d/2)**2
    return projected_area

def calculate_drag_force(reynolds_number,
    drag_coefficient,
    dynamic_viscosity,
    relative_velocity,
    air_density,
    projected_area,
    d):
    """Calculates drag force using the drag equation or stokes law depending on the calculated Reynolds number"""
    if reynolds_number < 1:
        drag_force = 3 * pi * dynamic_viscosity * relative_velocity * d
    else:
        drag_force = ( 0.5 * air_density * abs(relative_velocity) * relative_velocity * projected_area * drag_coefficient)
    return drag_force

def calculate_droplet_mass(droplet_volume, water_density):
    """calculates droplets mass in kg based on density and volume assuming the droplet is a perfect sphere"""
    droplet_mass = droplet_volume * water_density
    return droplet_mass

def calculate_buoyancy_force(air_density, droplet_volume):
    """calculates buoyancy force in N based on the weight of air displaced by the droplet,
     acts in the positive z-direction"""
    buoyancy_force = air_density * droplet_volume * 9.81
    return buoyancy_force

def calculate_force_of_gravity(droplet_mass):
    """calculates gravity force in N based on the mass of the droplet"""
    force_of_gravity = droplet_mass * 9.81
    return force_of_gravity

def calculate_acceleration_x_direction(drag_force_x_direction, mass):
    """calculates acceleration x direction in m/s^2"""
    acceleration_x_direction = drag_force_x_direction / mass
    return acceleration_x_direction


def calculate_acceleration_y_direction(drag_force_y_direction, mass):
    """calculates acceleration y direction in m/s^2"""
    acceleration_y_direction = drag_force_y_direction / mass
    return acceleration_y_direction

def calculate_acceleration_z_direction(drag_force_z_direction, mass, buoyancy_force, force_of_gravity):
    """calculates acceleration y direction in m/s^2"""
    sum_of_forces = drag_force_z_direction - force_of_gravity + buoyancy_force
    acceleration_z_direction = sum_of_forces / mass
    return acceleration_z_direction

def calculate_new_velocity(time_step, acceleration, current_velocity):
    """calculates new velocity in m/s"""
    new_velocity = time_step * acceleration + current_velocity
    return new_velocity

def calculate_new_position(time_step,velocity, current_position):
    """calculates new position in m"""
    new_position = time_step * velocity +  current_position
    return new_position


#  heat transfer calculations
# lumped droplet model
def calculate_droplet_area(d):
    """calculates droplet surface area in m^2"""
    droplet_area = d ** 2 * pi
    return droplet_area
def calculate_prandtl_number(dynamic_viscosity_of_air, specific_heat_of_air, k_air):
    """calculates prandtl number """
    return dynamic_viscosity_of_air * specific_heat_of_air / k_air

def calculate_nusselt_number(prandtl_number, reynolds_number ):
    """calculates Nusselt number using a curve fit equation"""
    nusselt_number = (2 + 0.66 / ((1 + (0.84 * prandtl_number ** (1 / 6)) ** 3) ** (1 / 3))
                      *(((reynolds_number * prandtl_number) ** 1.7) / (
                    1 + (prandtl_number * reynolds_number) ** 1.2)))
    return nusselt_number

def calculate_thermal_conductivity_of_water(T_surface):
    """calculates thermal conductivity of water in W/m*K with a curve fit equation"""
    k_water = 0.5579 +0.0021* T_surface-(T_surface**2)*10**(-5)
    return k_water

def calculate_thermal_conductivity_of_air(T_film):
    """calculates the thermal conductivity of air in W/m*K using a curve fit equation"""
    k_air = 0.024139 + T_film * 7.655*10**(-5)
    return k_air

def calculate_convective_heat_transfer_coefficient(nusselt_number, k_air, d):
    """calculates the convective heat transfer coefficient between the droplet and the air in W/m^2*K"""
    h =  nusselt_number * k_air / d
    return h

# layered shell heat transfer
def calculate_biot_number(h, d, k_water):
    """calculates the Biot number"""
    biot_number = h * (d/6) / k_water
    return biot_number

def calculate_convective_heat_transfer(h, air_temperature, layer_temperature, droplet_area, time_step):
    """calculates the outer layer heat transfer rate through convection in J"""
    heat_transfer_rate = (air_temperature - layer_temperature) * h * droplet_area
    return heat_transfer_rate* time_step

def calculate_conductive_heat_transfer(neighbor_temperature, current_temperature, k_water, area, time_step, dr):
    """returns the convective heat transfer coefficient"""
    heat_transfer_rate = k_water*(neighbor_temperature - current_temperature)* area/dr
    return heat_transfer_rate * time_step


# evaporation calculations
def calculate_diffusion_coefficient(T_film, pressure):
    """calculates the diffusion coefficient"""
    return 2.4e-5 * ((T_film + 273.15)/300) ** 2.072 * (101325/pressure)

def calculate_schmidt_number(dynamic_viscosity, air_density, diffusion_coefficient):
    """calculates the schmidt number"""
    return dynamic_viscosity / (air_density * diffusion_coefficient)

def calculate_sherwood_number(schmidt_number, reynolds_number):
    """calculates the sherwood number"""
    return schmidt_number ** (1/3) *0.6* reynolds_number **.5 + 2

def calculate_surface_saturation_vapour_pressure(
    droplet_temperature):
    """returns saturation vapour pressure at droplet surface (Pa)."""
    return 610.788 * e ** ((17.27 * droplet_temperature)/ (droplet_temperature + 237.3))

def calculate_ambient_vapour_pressure(air_temperature,relative_humidity):
    """returns ambient water vapour pressure (Pa)."""
    return relative_humidity* 610.788* e ** ((17.27 * air_temperature) / (air_temperature + 237.3))

def calculate_humidity_ratio(pressure, vapour_pressure):
    """return the relative humidity ratio"""
    return 0.622*(vapour_pressure/ (pressure - vapour_pressure))

def calculate_humid_air_density(specific_gas_constant_of_water_vapour,specific_gas_constant_of_air, air_pressure,vapour_partial_pressure, T_film ):
    """Returns the humid density of air in kg/m^3"""
    dry_air_partial_pressure = air_pressure - vapour_partial_pressure
    humid_air_density = dry_air_partial_pressure / (specific_gas_constant_of_air * (T_film +273.15)) + vapour_partial_pressure / (specific_gas_constant_of_water_vapour * (T_film + 273.15))
    return humid_air_density

def calculate_mass_transfer_coefficient(sherwood_number, diffusion_coefficient, diameter):
    """returns the mass transfer coefficient in m/s"""
    return sherwood_number * (diffusion_coefficient / diameter)

def calculate_mass_loss_per_second_humidity(k_m, diameter, surface_humidity_ratio, ambient_humidity_ratio,humid_air_density):
    """calculates the mass loss per second in kg/s"""
    return (pi * diameter**2 * k_m
            * (surface_humidity_ratio
            - ambient_humidity_ratio)
            * humid_air_density)

def calculate_mass_loss_per_second_heat(h, d, air_temperature, droplet_temperature):
    """calculates the mass loss per second in kg/s"""
    return h * pi * d**2 * (air_temperature - droplet_temperature) / 2.26e6

def calculate_mass_loss(time_step, mass_loss_per_second_humidity, mass_loss_per_second_heat):
    """returns mass loss in kg"""
    if abs(mass_loss_per_second_humidity) > abs(mass_loss_per_second_heat):
        return mass_loss_per_second_heat * time_step
    else:
        return mass_loss_per_second_humidity * time_step

def calculate_new_mass(droplet_mass, mass_loss ):
    """calculates the new mass in kg"""
    return droplet_mass - mass_loss

def update_diameter(new_mass, water_density):
    """updates the diameter of the droplet using the new calculated mass"""
    new_diameter = ((6 * new_mass) / ( pi * water_density))**(1/3)
    return new_diameter

def run_droplet_simulation(d: list, number_of_runs, time_step, number_of_shells):
    """returns a pandas data frame of simulated water droplet time steps given a list of diameters, number of runs, time step duration, number of shells """

    # constants
    specific_gas_constant_of_air = 287.05  # kj/kg*K
    specific_heat_of_water = 4180  # J/kg*k
    specific_heat_of_air = 1005 # J/kg*k
    specific_gas_constant_or_water_vapour = 461.5  # j/kg*K
    original_diameter = d #m

    x_velocity, y_velocity, z_velocity = look_up_initial_droplet_velocity()

    history = []

    # creating initial shell layer temperature
    shell_temperatures = []
    for j in range(number_of_shells):
        temperature = look_up_droplet_temperature()
        shell_temperatures.append(temperature)

    # creating initial droplet state which updates after every time step
    state = {
        'Time': 0,
        'x': 11.45,
        'y': 0,
        'z': 0,
        'x Velocity': x_velocity,
        'y Velocity': y_velocity,
        'z Velocity': z_velocity,
        'Center Droplet Temperature': look_up_droplet_temperature(),
        "Droplet Edge Temperature": look_up_droplet_temperature(),
        "Mean Droplet Temperature": look_up_droplet_temperature(),
        'Droplet Diameter': d
    }

    # main loop
    for step in range(number_of_runs-1):

        time = state['Time']
        time +=  time_step # updating time after every iteration


        # wind tunnel spatial conditions
        flow_x_velocity, flow_y_velocity, flow_z_velocity = look_up_flow_velocities(state['x'], state['y'], state['z'])
        pressure = look_up_pressure(state['x'], state['y'], state['z'])
        relative_humidity = look_up_relative_humidity(state['x'], state['y'], state['z'])
        air_temperature = look_up_air_temperature(state['x'], state['y'], state['z'])

        # droplet parameters and position
        x_velocity, y_velocity, z_velocity = state['x Velocity'], state['y Velocity'], state['z Velocity']
        x, y, z = state['x'], state['y'], state['z']
        T_center = state['Center Droplet Temperature']
        T_surface = state['Droplet Edge Temperature']
        T_mean = state["Mean Droplet Temperature"]
        d = state['Droplet Diameter']
        radius = d/2


        local_conditions = {
            'Pressure': pressure,
            'Humidity': relative_humidity,
            'Air Temperature': air_temperature,
            'Flow x Velocity': flow_x_velocity,
            'Flow y Velocity': flow_y_velocity,
            'Flow z Velocity': flow_z_velocity,
        }

        # calculations
        # force and kinematic
        v_rel_x = calculate_relative_velocity(flow_x_velocity, x_velocity)
        v_rel_y = calculate_relative_velocity(flow_y_velocity, y_velocity)
        v_rel_z = calculate_relative_velocity(flow_z_velocity, z_velocity)

        v_rel = calculate_relative_velocity_magnitude(v_rel_x, v_rel_y, v_rel_z)
        T_film = calculate_film_temperature(air_temperature, T_surface)

        surface_vapour_partial_pressure = calculate_surface_saturation_vapour_pressure(T_surface)
        ambient_vapour_partial_pressure = calculate_ambient_vapour_pressure(air_temperature, relative_humidity)
        ambient_humidity_ratio = calculate_humidity_ratio(pressure, ambient_vapour_partial_pressure)
        surface_humidity_ratio = calculate_humidity_ratio(pressure,surface_vapour_partial_pressure)
        humid_air_density = calculate_humid_air_density(specific_gas_constant_or_water_vapour,specific_gas_constant_of_air,pressure, ambient_vapour_partial_pressure, T_film)

        droplet_volume = calculate_droplet_volume(d)
        air_density = calculate_air_density(pressure, specific_gas_constant_of_air, air_temperature)
        water_density = calculate_water_density(T_mean)
        droplet_mass = calculate_droplet_mass(droplet_volume, water_density)
        droplet_area = calculate_droplet_area(d)
        droplet_projected_area = calculate_projected_area(d)


        k_air = calculate_thermal_conductivity_of_air(T_film)
        k_water = calculate_thermal_conductivity_of_water(T_mean)

        dynamic_viscosity_film_temp = calculate_dynamic_viscosity(T_film)
        dynamic_viscosity_of_air = calculate_dynamic_viscosity(air_temperature)
        reynolds_number = calculate_reynolds_number(humid_air_density, v_rel, d, dynamic_viscosity_film_temp)
        drag_coefficient = calculate_drag_coefficient(reynolds_number)

        drag_force_x_direction = calculate_drag_force(reynolds_number,drag_coefficient,dynamic_viscosity_film_temp,v_rel_x, humid_air_density, droplet_projected_area, d)
        drag_force_y_direction = calculate_drag_force(reynolds_number, drag_coefficient, dynamic_viscosity_film_temp,v_rel_y, humid_air_density, droplet_projected_area, d)
        drag_force_z_direction = calculate_drag_force(reynolds_number, drag_coefficient, dynamic_viscosity_film_temp, v_rel_z, humid_air_density, droplet_projected_area, d)
        force_of_gravity = calculate_force_of_gravity(droplet_mass)
        buoyancy_force = calculate_buoyancy_force(humid_air_density, droplet_volume)
        x_acceleration = calculate_acceleration_x_direction(drag_force_x_direction, droplet_mass)
        y_acceleration = calculate_acceleration_y_direction(drag_force_y_direction, droplet_mass)
        z_acceleration = calculate_acceleration_z_direction(drag_force_z_direction, droplet_mass, buoyancy_force, force_of_gravity)
        d_ratio = d/original_diameter

        new_x_velocity = calculate_new_velocity(time_step, x_acceleration, x_velocity)
        new_y_velocity = calculate_new_velocity(time_step, y_acceleration, y_velocity)
        new_z_velocity = calculate_new_velocity(time_step, z_acceleration, z_velocity)

        new_x = calculate_new_position(time_step, x_velocity,x)
        new_y = calculate_new_position(time_step, y_velocity, y)
        new_z = calculate_new_position(time_step, z_velocity, z)

        force_kinetic_calculations = {
            'Relative x Velocity': v_rel_x,
            'Relative y Velocity': v_rel_y,
            'Relative z Velocity': v_rel_z,
            'Relative Velocity Magnitude': v_rel,
            'Air Density': air_density,
            'Water Density': water_density,
            'Droplet Mass': droplet_mass,
            'Droplet Area': droplet_area,
            'Droplet Projected Area': droplet_projected_area,
            'Film Temperature': T_film,
            'Thermal Conductivity of Air': k_air,
            "Thermal Conductivity of Water": k_water,
            "Dynamic Viscosity Film Temp": dynamic_viscosity_film_temp,
            "Dynamic Viscosity of Air": dynamic_viscosity_of_air,
            'Reynolds Number': reynolds_number,
            'Drag Coefficient': drag_coefficient,
            "Drag Force X Direction": drag_force_x_direction,
            "Drag Force Y Direction": drag_force_y_direction,
            'Drag Force Z Direction': drag_force_z_direction,
            "Force of gravity": force_of_gravity,
            'Buoyancy Force': buoyancy_force,
            'X Acceleration': x_acceleration,
            'Y Acceleration': y_acceleration,
            'Z Acceleration': z_acceleration,
            'Droplet Diameter Ratio': d_ratio

        }
        # heat transfer

        prandtl_number = calculate_prandtl_number(dynamic_viscosity_film_temp, specific_heat_of_air, k_air)
        nusselt_number = calculate_nusselt_number(prandtl_number,reynolds_number)

        # heat transfer shell module
        dr = radius/number_of_shells

        h = calculate_convective_heat_transfer_coefficient(nusselt_number, k_air, d)

        if step == 0:
            biot_number = calculate_biot_number(h, d, k_water)
        else:
            pass
        updated_biot_number = calculate_biot_number(h, d,k_water)
        new_shell_masses = []
        new_shell_temperatures = []
        for j in range(number_of_shells):

            r_center = (j + 0.5) * dr # not used but shows where the center of each layer is
            r_in = j * dr
            r_out = (j + 1) * dr

            volume = 4 / 3 * np.pi * (r_out ** 3 - r_in ** 3)
            shell_mass = water_density * volume

            area_in = 4 * np.pi * r_in ** 2
            area_out = 4 * np.pi * r_out ** 2

            if j == number_of_shells - 1:
                q_from_outer_layer = calculate_convective_heat_transfer(h, air_temperature,shell_temperatures[j],area_out, time_step )
                q_from_inner_layer = calculate_conductive_heat_transfer(shell_temperatures[j], shell_temperatures[j-1],k_water, area_in, time_step, dr)

            elif j == 0:
                q_from_outer_layer = calculate_conductive_heat_transfer(shell_temperatures[j+1], shell_temperatures[j], k_water, area_out, time_step, dr)
                q_from_inner_layer = 0
            else:
                q_from_outer_layer = calculate_conductive_heat_transfer(shell_temperatures[j+1],shell_temperatures[j],k_water, area_out, time_step, dr)
                q_from_inner_layer = calculate_conductive_heat_transfer(shell_temperatures[j], shell_temperatures[j-1], k_water, area_in,time_step, dr)


            total_heat_transfer = q_from_outer_layer - q_from_inner_layer

            new_shell_masses.append(shell_mass)
            new_shell_temperature = shell_temperatures[j] + (total_heat_transfer / (shell_mass * specific_heat_of_water))
            new_shell_temperatures.append(new_shell_temperature)


        shell_temperatures = np.array(new_shell_temperatures)
        shell_masses = np.array(new_shell_masses)

        # use updated biot number if you want heat transfer methods to change part way through and use biot number if you want the method to stay constant
        # If you want the biot number to change the method of heat transfer a better method to switch between the two methods need to be implemented in order to do so
        # right now the shell model is being calculated the whole time and when it switches it does not take into account the lumped models temperature
        # when switching, the shells need to be set to the current lumped model temperature

        if biot_number < 0.1: #uses the original biot number so the calculation method doesn't change part way through the simulation
            total_heat_transfer = calculate_convective_heat_transfer(h, air_temperature,T_surface, droplet_area, time_step)
            final_center_temperature = T_center + total_heat_transfer / (droplet_mass * specific_heat_of_water)
            new_T_center = final_center_temperature
            new_T_surface = final_center_temperature
            T_mean = final_center_temperature
        else:
            final_center_temperature = shell_temperatures[0]
            final_edge_temperature = shell_temperatures[-1]
            new_T_center =  final_center_temperature
            new_T_surface = final_edge_temperature
            T_mean = np.sum(shell_masses * shell_temperatures) / np.sum(shell_masses)

        heat_transfer_calculations = {
            "Prandtl Number": prandtl_number,
            "Nusselt Number": nusselt_number,
            "Convective Heat Transfer Coefficient": h,
            "Original Biot Number": biot_number,
            "Biot Number": updated_biot_number,
            "Center Shell Temperature": new_T_center,
            "Outer Shell Temperature": new_T_surface,
            "Mean Droplet Temperature": T_mean
        }


        # evaporation
        diffusion_coefficient = calculate_diffusion_coefficient(T_film, pressure)
        schmidt_number = calculate_schmidt_number(dynamic_viscosity_film_temp,humid_air_density, diffusion_coefficient)
        sherwood_number = calculate_sherwood_number(schmidt_number, reynolds_number)

        # surface_vapour_partial_pressure = calculate_surface_saturation_vapour_pressure(T_surface)
        # ambient_vapour_partial_pressure = calculate_ambient_vapour_pressure(air_temperature, relative_humidity)
        # ambient_humidity_ratio = calculate_humidity_ratio(pressure, ambient_vapour_partial_pressure)
        # surface_humidity_ratio = calculate_humidity_ratio(pressure, surface_vapour_partial_pressure)
        # humid_air_density = calculate_humid_air_density(specific_gas_constant_or_water_vapour,
        #                                                 specific_gas_constant_of_air, pressure,
        #                                                 ambient_vapour_partial_pressure, T_film)

        k_m = calculate_mass_transfer_coefficient(sherwood_number, diffusion_coefficient, d)
        mass_loss_per_second_humidity = calculate_mass_loss_per_second_humidity(k_m, d,surface_humidity_ratio, ambient_humidity_ratio, humid_air_density)
        mass_loss_per_second_heat = calculate_mass_loss_per_second_heat(h,d, air_temperature, T_surface)
        mass_loss = calculate_mass_loss(time_step, mass_loss_per_second_humidity, mass_loss_per_second_heat)

        new_mass =  calculate_new_mass(droplet_mass, mass_loss)
        new_mass = max(new_mass, 0)

        new_droplet_diameter = update_diameter(new_mass, water_density)

        mass_transfer_calculations = {
            "Diffusion Coefficient": diffusion_coefficient,
            "Schmidt Number": schmidt_number,
            "Sherwood Number": sherwood_number,
            "Saturation Vapour Pressure": surface_vapour_partial_pressure,
            "Ambient Vapour Pressure": ambient_vapour_partial_pressure,
            "Mass Transfer Coefficient": k_m,
            'Mass loss per second_humidity': mass_loss_per_second_humidity,
            'Mass loss per second_heat': mass_loss_per_second_heat,
            'Mass loss': mass_loss

        }

        history.append(
            {"Time": time_step,
                        **state,
                        **local_conditions,
                        **force_kinetic_calculations,
                        **heat_transfer_calculations,
                        **mass_transfer_calculations,
             }
        )

        state = {
    'Time': time,
    'x': new_x,
    'y': new_y,
    'z': new_z,
    'x Velocity': new_x_velocity,
    'y Velocity': new_y_velocity,
    'z Velocity': new_z_velocity,
    'Center Droplet Temperature': new_T_center,
    'Droplet Edge Temperature': new_T_surface,
    'Mean Droplet Temperature': T_mean,
    'Droplet Diameter': new_droplet_diameter
}

    return pd.DataFrame(history)


# plotting
def plot_position_vs_position(history_data_frame, x_axis, y_axis):
    """plots a 2D plot of one position vs another, where x,y and z are given the values 1, 2, and 3 respectively"""
    fig, ax = plt.subplots()
    position_map = {
        1: ('x','X Position'),
        2: ('y','Y Position'),
        3: ('z','Z Position'),
    }
    if x_axis in position_map:
        x_column_name, x_axis_name = position_map[x_axis]
    else:
        print('X axis not recognized \n Input 1, 2, or three for the X axis ')
    if y_axis in position_map:
        y_column_name, y_axis_name = position_map[y_axis]
    else:
        print('Y axis not recognized \n Input 1, 2, or three for the Y axis ')

    ax.plot(history_data_frame[x_column_name], history_data_frame[y_column_name], label=f"{x_axis_name} vs {y_axis_name}")
    x = history_data_frame[x_column_name]
    y = history_data_frame[y_column_name]

    x_min, x_max = x.min(), x.max()
    y_min, y_max = y.min(), y.max()

    range_max = max(
        x_max - x_min,
        y_max - y_min
    )

    x_mid = (x_max + x_min) / 2
    y_mid = (y_max + y_min) / 2

    ax.set_xlim(
        x_mid - range_max / 2,
        x_mid + range_max / 2
    )

    ax.set_ylim(
        y_mid - range_max / 2,
        y_mid + range_max / 2
    )

    ax.set_aspect('equal')

    ax.set_title(f"{x_axis_name} vs {y_axis_name}")
    ax.set_aspect('equal', adjustable='box')
    ax.set_ylabel(f"{y_axis_name} (m)")
    ax.set_xlabel(f"{x_axis_name} (m)")
    ax.grid(True)
    ax.legend()
    plt.show()

    return fig

def plot_vs_position(history_data_frames, labels, y_axis_name):
    """plots any calculated parameter vs x position"""
    fig, ax = plt.subplots(figsize=(12,6))
    for df, label in zip(history_data_frames, labels):
        ax.plot(df['x'],
                df[y_axis_name],
                label=label)

    ax.set_xlabel("Distance (m)", fontsize = 14)
    ax.set_ylabel(y_axis_name, fontsize = 14)
    ax.set_title(f"{y_axis_name} vs Horizontal Distance",fontsize = 16, fontweight='bold')
    ax.set_xlim(left=-0.75)
    # ax.axvline(11.45, color='red', linestyle='-', linewidth=2, label = 'Test Section')
    ax.grid(True)
    ax.legend()
    plt.show()

    return fig

def plot_vs_time(history_data_frames, labels, y_axis_name):
    """plots any parameter vs time"""

    fig, ax = plt.subplots(figsize=(12,6))
    for df, label in zip(history_data_frames, labels):

        ax.plot(df['Time'],
                df[y_axis_name],
                label=label)

    ax.set_xlabel("Time (s)", fontsize = 14, fontweight='bold')
    ax.set_ylabel(y_axis_name,fontsize = 14, fontweight='bold')
    ax.set_title(f"{y_axis_name} vs Time",fontsize = 16, fontweight='bold')
    ax.set_xlim(left=0)
    ax.grid(True)
    ax.legend()
    plt.show()

    return fig

if __name__ == '__main__':
    # diameters = [458.9, 522.8, 595.6, 678.6, 773.1, 880.8, 1003.5, 1143.2]
    diameters = [595.6]
    time_step = 0.00001
    number_of_runs = 30000
    number_of_shells = 15
    history_data_frames = []
    labels =[]
    for diameter in diameters:
        history_data_frame = run_droplet_simulation(diameter * 1e-6, number_of_runs, time_step, number_of_shells)
        print(history_data_frame)
        history_data_frames.append(history_data_frame)
        label = f'Diameter = {diameter}µm'
        labels.append(label)
    plot_vs_time(history_data_frames, labels, 'Mean Droplet Temperature')
    plot_vs_position(history_data_frames, labels, 'Mean Droplet Temperature')
    plot_vs_time(history_data_frames, labels, 'x Velocity')
    plot_vs_position(history_data_frames, labels, 'x Velocity')
    plot_vs_position(history_data_frames, labels, 'Relative Velocity Magnitude')
    plot_vs_position(history_data_frames, labels, 'Relative x Velocity')
    plot_vs_position(history_data_frames, labels, 'Droplet Diameter Ratio')
    # plot_vs_position(history_data_frames, labels, 'Biot Number')


