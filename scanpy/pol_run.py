import mcapic, thor, time
import numpy as py

#acqtime, filebase, step = 100, "100412/100412a", 10
#acqtime, filebase, step = 100, "100412/100412b", -10
# for c shift 3 inwards in horizontal positioner, add beta blocker
#acqtime, filebase, step = 100, "100412/100412c", -10
# for d remove beta blocker
#acqtime, filebase, step = 100, "100412/100412d", -10
# for e shift outwards 6 in horizontal positioner
#acqtime, filebase, step = 100, "100412/100412e", 10

# for 100413a shift inwards 2 in horizontal positioner
#acqtime, filebase, step = 100, "100413/100413a", 10
# for 100413b, no changes
#acqtime, filebase, step = 100, "100413/100413b", 10
# for 100413c shift inwards 2 in horizontal positioner @ 19
#acqtime, filebase, step = 100, "100413/100413c", -10
# for 100413d no change, run forwards
#acqtime, filebase, step = 100, "100413/100413d", 10

# for 100414a no change, run forwards
#acqtime, filebase, step = 100, "100414/100414a", 10
#acqtime, filebase, step = 100, "100414/100414b", 10
#acqtime, filebase, step = 100, "100414/100414c", -10
#acqtime, filebase, step = 100, "100414/100414d", -10
#acqtime, filebase, step = 100, "100414/100414e", -10

#acqtime, filebase, step = 30, "100415/100415a", 10
#acqtime, filebase, step = 60, "100416/100416a", 10
#acqtime, filebase, step = 60, "100416/100416b", 10
#acqtime, filebase, step = 40, "100416/100416c", 15
#acqtime, filebase, step = 40, "100419/100419a", -15
#acqtime, filebase, step = 20, "100419/100419b", 15
#acqtime, filebase, step = 20, "100419/100419c", -15
#acqtime, filebase, step = 40, "100419/100419d", 10
#acqtime, filebase, step = 40, "100419/100419e", -10
#acqtime, filebase, step = 40, "100419/100419f", 10
#acqtime, filebase, step = 40, "100419/100419g", 10
#acqtime, filebase, step = 40, "100419/100419h", 15
#acqtime, filebase, step = 40, "100419/100419h", 15
#acqtime, filebase, step = 20, "100420/100420a", 15
#acqtime, filebase, step = 40, "100420/100420b", -10
#acqtime, filebase, step = 40, "100420/100420c", 10
#acqtime, filebase, step = 40, "100420/100420d", -10
#acqtime, filebase, step = 40, "100420/100420e", -10
#acqtime, filebase, step = 40, "100420/100420f", -10
#acqtime, filebase, step = 20, "100421/100421a", 15
#acqtime, filebase, step = 20, "100421/100421b", -15
#acqtime, filebase, step = 30, "100421/100421c", -10
#acqtime, filebase, step = 30, "100421/100421d", 10
#acqtime, filebase, step = 30, "100421/100421e", -10
#acqtime, filebase, step = 30, "100421/100421f", 15
#acqtime, filebase, step = 20, "100421/100421g", 15
#acqtime, filebase, step = 10, "100421/100421h", 15
#acqtime, filebase, step = 10, "100430/100430a", 15
#acqtime, filebase, step = 10, "100430/100430b", -15
#acqtime, filebase, step = 10, "100430/100430c", 15
#acqtime, filebase, step = 10, "100430/100430d", 15
#acqtime, filebase, step = 10, "100430/100430e", -15
#acqtime, filebase, step = 10, "100430/100430f", -15
#acqtime, filebase, step = 10, "100430/100430g", 15

#acqtime, filebase, step = 10, "100503/100503a", 15
#acqtime, filebase, step = 10, "100503/100503b", -15
#acqtime, filebase, step = 10, "100503/100503c", 15
#acqtime, filebase, step = 20, "100503/100503d", 10
#acqtime, filebase, step = 20, "100503/100503e", -10
#acqtime, filebase, step = 20, "100503/100503f", -10
#acqtime, filebase, step = 20, "100503/100503g", 10
#acqtime, filebase, step = 20, "100503/100504a", -10

#acqtime, filebase, step = 25, "100609/100609a", -15

acqtime, filebase, step = 20, "test", -20

if (step > 0):
    angles = py.arange(0.0, 360+1, step)
elif (step < 0):
    angles = py.arange(360.0, -1, step)
else:
    print "Cannot have step = 0"
print "This run will take ", len(angles)*(20+acqtime*1.05+3+2)/60.0, " minutes"

# move to initial position
thor.move(angles[0])
time.sleep(20)

# loop over the requested rotation angles
for pos in angles:
    print "Taking data at angle = ", pos
    thor.move(pos)
    time.sleep(3)
    mcapic.acquire(acqtime)
    spos = "%03d" % pos
    filename = filebase+"_r"+spos+".pha"
    print "Saving with filename = ", filename
    keywords = [['ROTATION', pos, 'Angle of rotation stage'],
                ['XRAY_SRC', 'Varian VF-50J-Rh/S', 'X-ray source name'],
                ['SRC_SN', '71053-6W', 'X-ray tube serial number'],
                ['SRC_HV', 6.9, 'X-ray tube high voltage (kV)'],
                ['SRC_FILA', 2.51, 'X-ray tube filament current (A)'],
                ['SRC_CURR', 1.01, 'X-ray tube current (mA)'],
                ['POLARIZA', 'HOPG crystal grade ZYB', 'Polarization analyzer']]
#                ['POLARIZA', 'Ge (220) crystal', 'Polarization analyzer']]
    mcapic.save(filename, keywords=keywords)

print "Done."
