### Assignment 2 Report

#### Ex.1

I set the parallelogram located at plane z=0. Then I created a function "ray_parallelogram_intersection" to check if a ray intersects with an arbitrary parallelogram. The input is ray origin e, ray direction d, parallelogram origin a and two edges u and v. I use carmer's rule to compute the intersection, and judge if it's valid accroding to its value.

![image-20211007094221190](/Users/hehanlin/Library/Application Support/typora-user-images/image-20211007094221190.png)

Using the function above, I got the result of parallelogram in orthographic rays as follow.

<img src="/Users/hehanlin/Library/Application Support/typora-user-images/image-20211007094445700.png" alt="image-20211007094445700" style="zoom:50%;" />

I  also created a function to check if a ray intersects with an arbitrary sphere. The input is ray origin e, ray direction d, sphere center and sphere radius.

![image-20211007094849548](/Users/hehanlin/Library/Application Support/typora-user-images/image-20211007094849548.png)

In the perspective raytrace, if the parallelogram plane is parallel with the screen, then no matter how we move the perspective the shape won't change (only size will change). However for sphere, if we change the perspective, the shape we get will change, too. As we can see following, the parallelogram is still a parallelogram, but the projection of sphere change from circle to ellipse.

<img src="/Users/hehanlin/Library/Application Support/typora-user-images/image-20211007100052425.png" alt="image-20211007100052425" style="zoom: 67%;" />

Be careful, that doesn't mean for plane, the projection won't change shape in perspective raytrace. If the plane is not parallel with the screen, then the projection will change too, as follow.

<img src="/Users/hehanlin/Library/Application Support/typora-user-images/image-20211007100845202.png" alt="image-20211007100845202" style="zoom:67%;" />

#### Ex.2

In the exercise 2, I make use of the "ray_sphere_intersection" function defined above to compute the perspective raytrace of an arbitrary sphere. We can change the RGB components to get different colors.

![image-20211007104527426](/Users/hehanlin/Library/Application Support/typora-user-images/image-20211007104527426.png)

Ambient is an important parameter deside the light of the dark area(the shade under the light source). With higher Ambient parameter, the shade area become brighter.

![image-20211007105400712](/Users/hehanlin/Library/Application Support/typora-user-images/image-20211007105400712.png)

The diffuse coefficient decides the brightness of the area lighted by light source.

![image-20211007110034039](/Users/hehanlin/Library/Application Support/typora-user-images/image-20211007110034039.png)

The specular coefficient decides the brightness of the highlight area.

![image-20211007110526949](/Users/hehanlin/Library/Application Support/typora-user-images/image-20211007110526949.png)

And the phong exponent decides the area size of highlight spot.

![image-20211007111015056](/Users/hehanlin/Library/Application Support/typora-user-images/image-20211007111015056.png)

I also make a parallelogram version shading, the parameters' influence is just as same as sphere one so I won't paste too many pictures.

<img src="/Users/hehanlin/Library/Application Support/typora-user-images/image-20211008222310078.png" alt="image-20211008222310078" style="zoom:50%;" />

<img src="/Users/hehanlin/Library/Application Support/typora-user-images/image-20211008222608135.png" alt="image-20211008222608135" style="zoom:50%;" />

#### Conclusion

In this assignment, I learned the basic ideas of how to implement different raytrace, including orthographic one and perspective one. I got the basic idea of shading, and the influences of different parameters.