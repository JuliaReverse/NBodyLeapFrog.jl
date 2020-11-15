

datestr = "June 20, 1988 - June 20, 2028"

# ==============================================================================
# Plot the orbits of the planets. The sun gets a little covered up by the orbit
# of Mercury and Venus. 
# ==============================================================================

plt = plot([], [])
for i in 1:planets
    plot!(plt,a[i,:],b[i,:])
end
plt.title("Solar System, "+datestr,{'size':'14'});
plt.grid('on')
plt.axis('equal')
if (checkMaximallyPacked ==1)
    plt.legend(['Sun','Mercury','Venus','Earth','Mars','PlanetX','Jupiter','Saturn','Neptune','Uranus','Pluto']) 
else
    plt.legend(['Sun','Mercury','Venus','Earth','Mars','Jupiter','Saturn','Neptune','Uranus','Pluto'])     
end
# I don't plot Pluto because it's so far, but it's there
plt.xlabel("x (AU)", {'size':'14'});
plt.ylabel("y (AU)", {'size':'14'});
#savefig('NBodyOrbit10.png', bbox_inches='tight'); 
plt.show()

# ==============================================================================
#  Plot the change in total energy of the solar system. The energy oscillates due to machine error.
# If the time step is smaller, the change in total energy reduces. For 10 years, with a time step of .001 years for
# 10,000 iterations, the change in energy was order 10^7.
# ==============================================================================

plt.figure(1)
plt.plot(timeSpace,ChangeInE,'k')
# plt.legend(['Sun','Mercury','Venus','Earth','Mars','Jupiter','Saturn','Neptune','Uranus','Pluto'])
plt.title("Change in Energy, "+datestr,{'size':'14'});
plt.xlabel("dt (Year)", {'size':'14'});
plt.ylabel("$M_J * AU^2 / (2\pi*year)^2$", {'size':'14'});
plt.show();

# savefig('NBodyOrbit10Energy.png', bbox_inches='tight'); 

# ==============================================================================
#  Plot the change in momentum divided by total momentum of the solar system. Again, the error is due to machine error.
#  The machine error doesn't go away over time, but it oscillates, which I assume means the errors cancel out, most likely
#  due to the symmetry and oscillation of the orbits of the planets and that the leap frog method conserves energy. 
#  The change in momentum is not quite centered at zero, which suggests the barycenter is not exactly zero. 
# ==============================================================================

plt.figure(2)
plt.plot(timeSpace, hRatio,'k')
#plt.legend(['Sun','Mercury','Venus','Earth','Mars','Jupiter','Saturn','Neptune','Uranus','Pluto'])
plt.title("Angular Momentum, \n"+datestr,{'size':'14'});
plt.xlabel("dt (Year)", {'size':'14'});
plt.ylabel("Change in Momentum / Total Momentum", {'size':'14'});
plt.show();

#savefig('NBodyOrbit100Momentum.png', bbox_inches='tight'); 

