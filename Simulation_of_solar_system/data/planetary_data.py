'''
Planetary Data Library -- stored as Python dictionaries
'''

# gravitational constant
G_meters = 6.67430e-11       # m**3 / kg / s**2
G        = G_meters * 10**-9 # km**3/ kg / s**2
AU       =  149597870.7   # km

mercury = {
            'mass'  :   3.302e23,
            'x'     :  -2.21207e7,
            'y'     :  -6.68243e7,
            'u'     :   3.66622e1,
            'v'     :  -1.23026e1
            }

venus = {
		'name'            : 'Venus',
		# 'spice_name'      : 'VENUS BARYCENTER',
		# 'SPICE_ID'        : 2,
		'mass'  :   48.685e23,
        'x'     :  -1.08573e8,
        'y'     :  -3.78424e6,
        'u'     :   8.98465e-1,
        'v'     :  -3.51720e1,
		'mu'              : 3.2485859200000006E+05,
		'radius'          : 6051.8,
		"semi_maj_axis"   : 0.72333199 * AU, 
		"eccentricity"    : 0.00677323, 
		'traj_color'      : 'y'
		}

earth = {
		'name'            : 'Earth',
		# 'spice_name'      : 'EARTH',
		# 'SPICE_ID'        : 399,
		'mu'              : 5.972e24 * G,
		'radius'          : 6378.0,
		"semi_maj_axis"   : 1.00000011 * AU, 
		"eccentricity"    : 0.01671022, 
		# 'J2'              : 1.081874e-3,
		'traj_color'      : 'b',
        'mass'  :   5.97219e24,
        'x'     :  -2.62790e7,
        'y'     :   1.44510e8,
        'u'     :  -2.98305e1,
        'v'     :  -5.22046
		}

moon = {
		'name'            : 'Moon',
		'mass'            : 5.972e24,
		'mu'              : 5.972e24 * G,
		'radius'          : 1737.4,
		# 'J2'              : 1.081874e-3,
		"eccentricity"    : 0.0549, 
		'traj_color'      : 'b'
		}

mars = {
		'name'            : 'Mars',
		# 'spice_name'      : 'MARS BARYCENTER',
		# 'SPICE_ID'        : 4,
		'mass'  :   6.4171e23,
        'x'     :   2.06927e8,
        'y'     :  -3.56073e6,
        'u'     :   1.30431,
        'v'     :   2.62815e1,
		'mu'              : 4.282837362069909E+04,
		'radius'          : 3397.0,
		"semi_maj_axis"   : 1.52366231 * AU, 
		"eccentricity"    : 0.09341233, 
		'traj_color'      : 'r'
		}

jupiter = {
		'name'            : 'Jupiter',
		'mass'  :   1.8981e27,
        'x'     :   5.97841e8,
        'y'     :   4.38704e8,
        'u'     :  -7.89263,
        'v'     :   1.11503e1,
		'mu'              : 1.26686e8,
		'radius'          : 71490.0,   # km
		"semi_maj_axis"   : 5.20336301 * AU, 
		"eccentricity"    : 0.04839266, 
		'traj_color'      : 'C3'
}



saturn = {
	'name'            : 'Saturn',
	'mass'            : 568.34e24,
	'radius'          : 58232.0,
	'mu'              : 37.931e6,
	"semi_maj_axis"   : 9.53707032 * AU, 
	"eccentricity"    : 0.05415060, 
	'traj_color'      : 'C2'
}



sun = {
	'name'            : 'Sun',
	'mass'            : 1.989e30,
	'mu'              : 1.3271244004193938E+11,
	'radius'          : 695510.0,

}

bodies = [
	venus, earth, moon, mars, 
	jupiter, saturn, sun ]

for body in bodies:
	body[ 'diameter' ] = body[ 'radius' ] * 2
