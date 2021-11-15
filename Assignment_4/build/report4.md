### Assignment 4

#### Ex.1

I implemented triangle-ray intersection by solving the system $f(u,v)=p(t)$. As well as the ray_color function, render_scene function and all other related function used in assignment3. However, since we don't need reflection, refraction and depth of field in assignment4, I use a struct "config" to control the enablement of these features.

In the exercise one, I haven't implement a bah. It takes several seconds to produce "bunny", and seems impossible to render "dragon".

![image-20211027114041778](/Users/hehanlin/Library/Application Support/typora-user-images/image-20211027114041778.png)

<img src="/Users/hehanlin/Library/Containers/com.tencent.xinWeChat/Data/Library/Application Support/com.tencent.xinWeChat/2.0b4.0.9/e35fd6e067eccd4d95c9ab1f087f78f9/Message/MessageTemp/9e20f478899dc29eb19741386f9343c8/Image/2981635366552_.pic.jpg" alt="2981635366552_.pic" style="zoom:90%;" />

<div STYLE="page-break-after: always;"></div>

#### Ex.2

In this exercise, I implement a bvh tree by recursively split them in top-down order. And in each recursion, I order the nodes passed in along the longest axis. For intersect_box function, I implemented a linear solution to find if a ray intersect with the axis-aligned-box, shown as follow:

![image-20211027123658729](/Users/hehanlin/Library/Application Support/typora-user-images/image-20211027123658729.png)

After implementing the bah, I render "bunny" and "dragon" in seconds.

Time cost of "bunny", from 7.12s to 0.363s:

![image-20211027123936816](/Users/hehanlin/Library/Application Support/typora-user-images/image-20211027123936816.png)

Time cost of "dragno":

![image-20211027124048737](/Users/hehanlin/Library/Application Support/typora-user-images/image-20211027124048737.png)

<div STYLE="page-break-after: always;"></div>

Final result:

<img src="/Users/hehanlin/Library/Containers/com.tencent.xinWeChat/Data/Library/Application Support/com.tencent.xinWeChat/2.0b4.0.9/e35fd6e067eccd4d95c9ab1f087f78f9/Message/MessageTemp/9e20f478899dc29eb19741386f9343c8/Image/2991635366603_.pic.jpg" alt="2991635366603_.pic" style="zoom:90%;" />

If we enable reflection, refraction, and depth of filed, then we get:

<img src="/Users/hehanlin/Library/Containers/com.tencent.xinWeChat/Data/Library/Application Support/com.tencent.xinWeChat/2.0b4.0.9/e35fd6e067eccd4d95c9ab1f087f78f9/Message/MessageTemp/9e20f478899dc29eb19741386f9343c8/Image/3021635367180_.pic_hd.jpg" alt="3021635367180_.pic_hd" style="zoom:90%;" />

<div STYLE="page-break-after: always;"></div>

Plus shadow, all features enabled. It takes about 40s:

![image-20211027164137612](/Users/hehanlin/Library/Application Support/typora-user-images/image-20211027164137612.png)

<img src="/Users/hehanlin/Library/Containers/com.tencent.xinWeChat/Data/Library/Application Support/com.tencent.xinWeChat/2.0b4.0.9/e35fd6e067eccd4d95c9ab1f087f78f9/Message/MessageTemp/9e20f478899dc29eb19741386f9343c8/Image/3001635366984_.pic_hd.jpg" alt="3001635366984_.pic_hd" style="zoom:90%;" />

