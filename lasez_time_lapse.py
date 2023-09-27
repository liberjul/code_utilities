import time
import serial
from pywinauto import application
from pywinauto.keyboard import send_keys
from nanpy import SerialManager, ArduinoApi


interval = 15 # time interval in seconds
lamp_delay_before = 2
lamp_delay_after = 5
loop_time = 0.01 # time in seconds for time-checking loop. Smaller is more precise, but more energetically expensive.
wait_prop = 0.95 # proportion of interval to wait with time.sleep.
connection = SerialManager('COM3')
a = ArduinoApi(connection=connection)
a.pinMode(13, a.OUTPUT)
a.digitalWrite(13, a.HIGH)
start = time.time()
t1 = start
cycle = 1
while True:
    a.digitalWrite(13, a.LOW)
    time.sleep(lamp_delay_before*wait_prop)
    while (time.time() - start) < (interval*cycle):
        time.sleep(loop_time)
    send_keys("{F3}")
    t2 = t1
    t1 = time.time()
    print(F"Cycle time: {t1-t2}; Cycle number: {cycle}")
    send_keys("{ENTER}")
    a.digitalWrite(13, a.HIGH)
    cycle += 1
    time.sleep(lamp_delay_after*wait_prop)
    while (time.time() - start) < (interval*cycle - lamp_delay_before):
        time.sleep(loop_time)
