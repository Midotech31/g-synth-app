from rest_framework import generics, permissions, status
from rest_framework.response import Response
from rest_framework.throttling import ScopedRateThrottle
from rest_framework.views import APIView
from rest_framework_simplejwt.views import TokenObtainPairView

from apps.accounts.serializers import (
    ChangePasswordSerializer,
    LogoutSerializer,
    RegisterSerializer,
    UserSerializer,
)


class RegisterView(generics.CreateAPIView):


    permission_classes = (permissions.AllowAny,)
    serializer_class = RegisterSerializer
    throttle_classes = (ScopedRateThrottle,)
    throttle_scope = "register"

    def create(self, request, *args, **kwargs):
        serializer = self.get_serializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        user = serializer.save()
        return Response(UserSerializer(user).data, status=status.HTTP_201_CREATED)


class LoginView(TokenObtainPairView):


    throttle_classes = (ScopedRateThrottle,)
    throttle_scope = "login"


class LogoutView(APIView):


    def post(self, request):
        serializer = LogoutSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        serializer.save()
        return Response(status=status.HTTP_205_RESET_CONTENT)


class LogoutAllView(APIView):


    def post(self, request):
        revoked = request.user.revoke_all_tokens()
        return Response({"detail": "All sessions revoked.", "sessions_revoked": revoked})


class MeView(APIView):


    def get(self, request):
        return Response(UserSerializer(request.user).data)


class ChangePasswordView(APIView):


    def post(self, request):
        serializer = ChangePasswordSerializer(
            data=request.data, context={"request": request}
        )
        serializer.is_valid(raise_exception=True)
        serializer.save()
        return Response({
            "detail": "Password updated. All existing sessions were signed out.",
        })
