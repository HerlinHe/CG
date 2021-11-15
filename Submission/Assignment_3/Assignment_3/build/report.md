### Assignment 3

I implemented perspective camera, specular highlights, shadow, reflections, refractions and depth of field in this homework.

#### Ex.1

I stretch the pixel in this way, so no matter how we change the focal_length, the field_of_view is fixed and we can see the exactly same scene(don't consider the different caused by depth of field). And scale_x is proportionable with scale_y.

<img src="/Users/hehanlin/Library/Application Support/typora-user-images/image-20211019193905949.png" alt="image-20211019193905949" style="zoom: 67%;" />

#### Ex.2

I implemented shadow rays, and set offset so the shadow ray won't be blocked by itself.

<img src="/Users/hehanlin/Library/Application Support/typora-user-images/image-20211019194318914.png" alt="image-20211019194318914" style="zoom:67%;" />

<img src="/Users/hehanlin/Library/Application Support/typora-user-images/image-20211019194623932.png" alt="image-20211019194623932" style="zoom:40%;" />

If there is no offset, then many pixels may be dark. Just like following:

<img src="/Users/hehanlin/Library/Application Support/typora-user-images/image-20211019194525841.png" alt="image-20211019194525841" style="zoom:40%;" />

<div STYLE="page-break-after: always;"></div>

#### Ex.3

Implemented refract function to compute the refraction ray direction.

![image-20211019195206030](/Users/hehanlin/Library/Application Support/typora-user-images/image-20211019195206030.png)

Using recursion to compute the reflection color and refraction color.

![image-20211019195251810](/Users/hehanlin/Library/Application Support/typora-user-images/image-20211019195251810.png)

<div STYLE="page-break-after: always;"></div>

#### Ex.4

Randomly sampled ten rays from camera, and averaged their result to get pixel color. What I get is identical with the professor's sample.

<img src="/Users/hehanlin/Library/Application Support/typora-user-images/image-20211019195559242.png" alt="image-20211019195559242" style="zoom:50%;" />

If I change the focal length from 2.0 to 6.0, the rightmost sphere is much more sharp.

![image-20211019195958081](/Users/hehanlin/Library/Application Support/typora-user-images/image-20211019195958081.png)

There is a scene in wich more object are pictured. The rightmost sphere is much more clear than others. And refraction is working.

<img src="/Users/hehanlin/Library/Application Support/typora-user-images/image-20211019200119009.png" alt="image-20211019200119009" style="zoom:50%;" />

<div STYLE="page-break-after: always;"></div>

#### Ex.5

I implemented animation, you can check the gif in my submission.

![image-20211019214106336](/Users/hehanlin/Library/Application Support/typora-user-images/image-20211019214106336.png)

Uncomment the animation function, you can reproduce the result. It takes some time because there are recursions in both reflection and refraction. Or you can set the refraction color to zero so it will be much faster but without refraction pattern.