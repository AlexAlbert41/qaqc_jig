BAC = "Milano"

if BAC == "CERN":
    from utilities import TXP3510P

if BAC == "Milano":
    from utilities import keithley2231A

class PSU:
    def __init__(self):
        if BAC == "Milano":
            self.psu = keithley2231A.Keithley2231A()
        if BAC == "CERN":
            self.psu1 = TXP3510P('/dev/TTi-1')
            self.psu2 = TXP3510P('/dev/TTi-2')

    def getState(self):
        if BAC == "Milano":
            return self.psu.check_state()
        if BAC == "CERN":
            return self.psu1.getState() * self.psu2.getState()

    def powerOn(self):
        if BAC == "Milano":
            self.psu = keithley2231A.Keithley2231A()
            self.psu.set_state(1)
        if BAC == "CERN":
            self.psu1 = TXP3510P('/dev/TTi-1')            
            self.psu1.setCurrent(2.)
            self.psu1.setVoltage(5.3)
            self.psu1.powerOn()
            self.psu2 = TXP3510P('/dev/TTi-2')
            self.psu2.setCurrent(1.,1)
            self.psu2.setVoltage(12.,1)
            self.psu2.powerOn(1)
            self.psu2.setCurrent(1.,2)
            self.psu2.setVoltage(12.,2)
            self.psu2.powerOn(2)
            
    def powerOff(self):
        if BAC == "Milano":
            self.psu = keithley2231A.Keithley2231A()
            self.psu.set_state(0)
        if BAC == "CERN":
            self.psu1 = TXP3510P('/dev/TTi-1')            
            self.psu1.powerOff()
            self.psu2 = TXP3510P('/dev/TTi-2')
            self.psu2.powerOff(1)
            self.psu2.powerOff(2)
            
