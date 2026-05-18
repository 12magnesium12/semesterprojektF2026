import asyncio

from asyncua import Client, ua
import time
from datetime import datetime
import math 

url = "opc.tcp://192.168.8.101:4841"
namespace = "urn:B&R/pv/"


async def main():
    print(f"Connecting to {url} ...")
    async with Client(url=url,  watchdog_intervall=5) as client:
        # Find the namespace index
        nsidx = await client.get_namespace_index(namespace)
        print(f"Namespace Index for '{namespace}': {nsidx}")

        f = open("dataPID.csv", "w")
        f.write("t(s);x(m);u(N);theta(rad)\n")

        # Get the variable node for read / write
        theta_node = await client.nodes.root.get_child(f"0:Objects/4:PLC/6:Modules/6:::/6:Program/6:theta")
        x_node = await client.nodes.root.get_child(f"0:Objects/4:PLC/6:Modules/6:::/6:Program/6:x")
        u_node = await client.nodes.root.get_child(f"0:Objects/4:PLC/6:Modules/6:::/6:Program/6:u")
        K_scale_node = await client.nodes.root.get_child(f"0:Objects/4:PLC/6:Modules/6:::/6:Program/6:K_scale")


        is_recording = False
        start_time = 0.0

        while(True):
            K_scale = await K_scale_node.read_value()
            if(K_scale == 1):
                if not is_recording:
                    is_recording = True
                    start_time = time.perf_counter()  # Capture the exact start time
                    print("K_scale is 1. Recording started...")

                theta = await theta_node.read_value()
                x = await x_node.read_value()
                u = await u_node.read_value()

                current_time = time.perf_counter()
                t = current_time - start_time
                t_rounded = round(t, 4)
                print(t_rounded)

                theta_str = str(theta-math.pi).replace('.', ',')
                x_str = str(x).replace('.', ',')
                u_str = str(u).replace('.', ',')
                t_str = str(t_rounded).replace('.', ',')
                
                f.write(f"{t_str};{x_str};{u_str};{theta_str}\n")

                await asyncio.sleep(0.001)
                




        ## Calling a method
        #res = await client.nodes.objects.call_method(f"{nsidx}:ServerMethod", 5)
        #print(f"Calling ServerMethod returned {res}")


if __name__ == "__main__":
    asyncio.run(main())
