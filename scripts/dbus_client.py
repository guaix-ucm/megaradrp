import dbus



bus = dbus.SessionBus()

megara = bus.get_object('es.ucm.MegaraService', '/es/ucm/MegaraService/Instrument/MEGARA')
somecall = megara.expose("a", dbus_interface='es.ucm.InstrumentInterface')


megara_dev_iface = dbus.Interface(megara, dbus_interface='es.ucm.InstrumentInterface')
somecall = megara_dev_iface.expose("a")


conf = {'a': 'a', 'b': 2}
res = megara_dev_iface.configure(conf)

print res

for i in res:
    print i, res[i]

