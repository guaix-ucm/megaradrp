from gi.repository import GLib

import dbus
import dbus.service
import dbus.mainloop.glib


from simulation import create_instrument

class DemoException(dbus.DBusException):
    _dbus_error_name = 'com.example.DemoException'


class MegaraObject(dbus.service.Object):
    def __init__(self, session_bus, path):
        super(MegaraObject, self).__init__(session_bus, path)
        self.instrument = create_instrument()

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


if __name__ == '__main__':
    import logging

    logging.basicConfig(level=logging.DEBUG)

    dbus.mainloop.glib.DBusGMainLoop(set_as_default=True)

    session_bus = dbus.SessionBus()
    name = dbus.service.BusName("es.ucm.MegaraService", session_bus)
    object = MegaraObject(session_bus, '/MegaraObject')

    mainloop = GLib.MainLoop()
    print "Running example service."

    mainloop.run()
