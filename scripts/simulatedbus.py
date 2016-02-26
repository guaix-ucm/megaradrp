from gi.repository import GLib

import dbus
import dbus.service
import dbus.mainloop.glib


from simulation import create_instrument
from megaradrp.simulation.factory import MegaraImageFactory

class DemoException(dbus.DBusException):
    _dbus_error_name = 'es.ucm.DemoException'


class MegaraInstrument(dbus.service.Object):
    def __init__(self, session_bus, path):
        super(MegaraInstrument, self).__init__(session_bus, path)
        self.instrument = create_instrument()


    @dbus.service.signal("es.ucm.MegaraInterface")
    def signal1(self, messages):
        pass

    @dbus.service.method("es.ucm.MegaraInterface")
    def emit_signal1(self):
        self.signal1('signal1')
        return 'signal1 emitted'

    @dbus.service.method("es.ucm.InstrumentInterface",
                         in_signature='a{sv}', out_signature='a{sv}')
    def configure(self, mode):
        return mode

    @dbus.service.method("es.ucm.InstrumentInterface",
                         in_signature='s', out_signature='s')
    def expose(self, name):
        return name

    @dbus.service.method("es.ucm.MegaraInterface",
                         in_signature='s', out_signature='as')
    def set_mode(self, mode):
        self.instrument.set_mode(mode)
        return [mode, "something",
                session_bus.get_unique_name()]

    @dbus.service.method("es.ucm.MegaraInterface",
                         in_signature='s', out_signature='as')
    def set_vph(self, name):
        self.instrument.wheel.select(name)

        return [name, "something",
                session_bus.get_unique_name()]

    @dbus.service.method("es.ucm.MegaraInterface",
                         in_signature='', out_signature='')
    def RaiseException(self):
        raise DemoException('The RaiseException method does what you might '
                            'expect')

    @dbus.service.method("es.ucm.MegaraInterface",
                         in_signature='s', out_signature='(ss)')
    def set_cover(self, mode):
        self.instrument.set_cover(mode)
        return ("Hello Tuple", " from example-service.py")

    @dbus.service.method("es.ucm.MegaraInterface",
                         in_signature='', out_signature='')
    def exit(self):
        mainloop.quit()


class OtherInstrument(dbus.service.Object):
    def __init__(self, session_bus, path):
        super(OtherInstrument, self).__init__(session_bus, path)
        self.instrument = create_instrument()

    @dbus.service.method("es.ucm.InstrumentInterface",
                         in_signature='a{sv}', out_signature='a{sv}')
    def configure(self, mode):
        return mode

    @dbus.service.method("es.ucm.InstrumentInterface",
                         in_signature='s', out_signature='as')
    def expose(self, name):


        return [name, "something",
                session_bus.get_unique_name()]


class FactoryObject(dbus.service.Object):
    def __init__(self, session_bus, path):
        super(FactoryObject, self).__init__(session_bus, path)

        self.factory = MegaraImageFactory()

    @dbus.service.method("es.ucm.MegaraInterface",
                         in_signature='', out_signature='')
    def exit(self):
        mainloop.quit()


if __name__ == '__main__':
    import logging

    logging.basicConfig(level=logging.DEBUG)

    dbus.mainloop.glib.DBusGMainLoop(set_as_default=True)

    session_bus = dbus.SessionBus()
    name = dbus.service.BusName("es.ucm.MegaraService", session_bus)
    object11 = MegaraInstrument(session_bus, '/es/ucm/MegaraService/Instrument/MEGARA')
    object12 = OtherInstrument(session_bus, '/es/ucm/MegaraService/Instrument/Other')
    object2 = FactoryObject(session_bus, '/es/ucm/MegaraService/Factory')

    mainloop = GLib.MainLoop()
    print "Running services..."

    mainloop.run()
