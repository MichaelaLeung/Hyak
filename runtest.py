def runtest():
    print('hello')
    currentDT = datetime.datetime.now()
    name = "smart_output_"+str(currentDT.month)+"/"+str(currentDT.day)+"-"+str(currentDT.hour)+":"+str(currentDT.minute)

    HERE = os.path.dirname(os.path.abspath(__file__))
    place = os.path.join(HERE, name)
    os.mkdir(place)
    print(place)
