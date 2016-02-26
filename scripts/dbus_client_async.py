import dbus

def handle_reply(r):

    print r
    print 'done'
    #mainloop.quit()

def handle_error(error):

    print 'some error'
    print error
    mainloop.quit()

def makecalls(bus):
    megara = bus.get_object('es.ucm.MegaraService', '/es/ucm/MegaraService/Instrument/MEGARA')
    # This will timeout
    somecall = megara.expose("a",
                             dbus_interface='es.ucm.InstrumentInterface',
                             reply_handler=handle_reply,
                             error_handler=handle_error
                             )
    print somecall


def handle_signal(message):
    print 'received mssg', message
    pass


def connect(bus):
    megara = bus.get_object('es.ucm.MegaraService', '/es/ucm/MegaraService/Instrument/MEGARA')
    # This will timeout
    megara.connect_to_signal('signal1', handle_signal, dbus_interface='es.ucm.MegaraInterface',)

def catchall_signal_handler(*args, **kwargs):
    print "Caught signal (in catchall handler) "
    for key in kwargs:
        print 'key=',key, 'val=',kwargs[key]

    for arg in args:
        print 'arg', arg


def catchall_signal1_signals_handler(messg):
    print 'received signal1 with msg', messg

def catchall_testservice_interface_handler(hello_string, dbus_message):
    print "es.ucm.MegaraInterface interface says " + hello_string
    print " when it sent signal " + dbus_message.get_member()

if __name__ == '__main__':

    from gi.repository import GLib
    import dbus
    import dbus.service
    import dbus.mainloop.glib

    import logging

    logging.basicConfig(level=logging.DEBUG)

    dbus.mainloop.glib.DBusGMainLoop(set_as_default=True)

    session_bus = dbus.SessionBus()

    makecalls(session_bus)

    connect(session_bus)

    session_bus.add_signal_receiver(catchall_signal_handler,
                                    interface_keyword='dbus_interface',
                                    member_keyword='member')

    session_bus.add_signal_receiver(catchall_signal1_signals_handler,
                                    dbus_interface = "es.ucm.MegaraInterface",
                                    signal_name = "signal1")

    session_bus.add_signal_receiver(catchall_testservice_interface_handler,
                                    dbus_interface = "es.ucm.MegaraInterface",
                                    message_keyword='dbus_message')

    mainloop = GLib.MainLoop()
    print "Running async example."

    mainloop.run()