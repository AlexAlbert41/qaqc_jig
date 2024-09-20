import pyvisa
import serial



##########################
class Keithley2231A():
    """Instrument class for Keithley 2231A

    Args:
        * portname (str): port name
        * channel (str): channel name

    """

    def __init__(self, portname='ASRL/dev/keithley2231A::INSTR', chName="CH1"):
        self.instr = pyvisa.ResourceManager('@py').open_resource(portname)
        self.chName = chName
        self.instr.write("SYSTem:REMote")
        self.instr.write("INST:SEL "+self.chName)

    def query(self, query):
        """Pass a query to the power supply"""
        print(self.instr.query(query).strip())

    def meas_V(self):
        """read set voltage"""
        volt = self.instr.query("MEAS:VOLT?").strip()
        return(float(volt))

    def meas_I(self):
        """measure current"""
        curr = self.instr.query("MEAS:CURR?").strip()
        return(float(curr))

    def set_V(self, value):
        """set voltage"""
        return(self.instr.write("APPL "+self.chName+","+str(value)))

    def set_state(self, value):
        """Set the PS state (0: OFF, 1: RUNNING)"""
        return(self.instr.write("OUTP "+str(value)))
    
    def check_state(self):
        """Check the PS state (0: OFF, 1: RUNNING)"""
        return(int(self.instr.query("OUTP?").strip()))
    
    def getChannel(self):
        return(self.chName)
    
    def selectChannel(self,chName):
        self.chName = chName
        self.instr.write("INST:SEL "+self.chName)
