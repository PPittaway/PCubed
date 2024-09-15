import serial
import time

class switchingValve(object):

    def connect(self, COM):
        self.valve = serial.Serial(port = COM, baudrate = 9600, timeout = .1)
        print("Switching valve connected on port " + COM)
    
    def switch(self):
        command = '<switch>'
        self.valve.write(bytes(command, 'utf-8'))
        time.sleep(2)
